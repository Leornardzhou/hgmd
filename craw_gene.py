#!/usr/bin/env python
# -*- coding=utf-8 -*-
import re
import os
import sys
import time
import random
import datetime
import argparse
import functools

import bs4
import requests

from proxy import Proxy

reload(sys)
sys.setdefaultencoding('utf8')

# Basic
LOGIN_URL = 'http://www.hgmd.cf.ac.uk/ac/validate.php'
GENE_TAG_URL = 'http://hgmdtrial.biobase-international.com/hgmd/pro/geneTag.php'
GENE_TOTAL_URL = 'http://hgmdtrial.biobase-international.com/hgmd/pro/gene.php?gene='
GENE_DETAIL_URL = 'http://hgmdtrial.biobase-international.com/hgmd/pro/all.php'
MUTATION_URL = 'http://hgmdtrial.biobase-international.com/hgmd/pro/mut.php?acc='

HEADERS = {
    'User-Agent':
    'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/65.0.3325.181 Safari/537.36'
}
EMAIL = '2120141047@mail.nankai.edu.cn'
PASSWORD = 'HGMD992140'


def now(time_fmt='%Y-%m-%d %H:%M:%S'):

    return datetime.datetime.now().strftime(time_fmt)


# 自定义的装饰器
def try_again(N=10):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            n = 0
            while n < N:
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    n += 1
                    print '[warn] `{}` failed {}st time: {} [{}]'.format(func.__name__, n, e, now())
                    time.sleep(random.randint(5, 10))
            if n == N:
                print '[error] `{}` failed too many times({}) [{}]'.format(func.__name__, N, now())
        return wrapper
    return decorator


class CrawHGMD(object):

    def __init__(self, args, proxy):
        self.args = args
        self.proxy = proxy

    def start(self):

        # # test here
        # print self.get_reports('CG126032')
        # print self.get_reports('CG1616105')
        # print self.get_reports('CD011580')
        # print self.get_reports('CS164742')
        # print self.get_reports('CX003357')  # meeting abstract
        # exit()

        with open(self.args['out'] + '.database.hg19.txt', 'w') as out_hg19, \
                open(self.args['out'] + '.database.hg38.txt', 'w') as out_hg38, \
                open(self.args['out'] + '.failed.genelist', 'w') as out_failed:

            title = '#HGMD_ID Gene_Symbol Chrom Start End Mutation Phenotype Reports PMIDs Class'.split()
            out_hg19.write('\t'.join(title) + '\n')
            out_hg38.write('\t'.join(title) + '\n')

            for gene in self.genelist:

                print '\033[31m> get_gene_detail start: {} at {}\033[0m'.format(gene, now())
                
                context = self.get_gene_detail(gene)

                if context is None or context == 'no':
                    status = '[warn] no_such_gene' if context == 'no' else 'failed'
                    out_failed.write('{}\t{}\n'.format(gene, status))
                    continue

                for each in context:
                    for ref in ('ref_hg19', 'ref_hg38'):
                        temp = re.split(r'Chr|:|-', each[ref])
                        if len(temp) == 3:
                            chrom = temp[1]
                            start = end = temp[2]
                        elif len(temp) == 4:
                            chrom, start, end = temp[1:4]
                        line = '{hgmd_id}\t{gene}\t{chrom}\t{start}\t{end}\t{mutation}\t{phenotype}\t{reports}\t{pmids}\t{display}\n'.format(
                            gene=gene, chrom=chrom, start=start, end=end, **each)
                        if ref == 'ref_hg19':
                            out_hg19.write(line)
                        elif ref == 'ref_hg38':
                            out_hg38.write(line)
                print '\033[31m> get_gene_detail done: {} at {}\033[0m'.format(gene, now())
                

    @try_again()
    def get_gene_detail(self, gene):

        contexts = []

        # 每个基因重新登录
        session, proxies = self.login()

        soup = self.get_soup(GENE_TOTAL_URL + gene, session, proxies)

        if 'no entries associated' in soup.select('table .center')[0].text:
            print '[warn] no such gene: {}'.format(gene)
            return 'no'

        counter = 0

        for each in soup.select('table.gene')[2].select('tr')[1:]:
            database = each.select('input[name="database"]')[0].attrs['value']
            print '> use database: {}'.format(database)
            data = {
                'gene': gene,
                'database': database
            }

            soup = self.get_soup(GENE_DETAIL_URL, session, proxies, method='POST', data=data)
            forms = soup.select('form[action="mut.php"]')
            for form in forms:

                # 每爬取20个条目后重新登录
                counter += 1
                if counter % 20 == 0:
                    session, proxies = self.login()

                hgmd_id = form.select('input')[-1].attrs['value']
                print '>> {} {} crawling HGMD: {}'.format(counter, gene, hgmd_id)

                tds = form.findParent('tr').select('td')
                mutation = '|'.join(each.text for each in tds[1:-4])
                display = tds[-4].text
                phenotype = tds[-3].text or '.'
                report_list, pmid_list = self.get_reports(hgmd_id, session, proxies)
                reports = '|'.join(report_list) or '.'
                pmids = '|'.join(pmid_list) or '.'

                ref_hg19 = ref_hg38 = 'Chr.:.-.'
                for span in tds[-1].select('span.gen'):
                    if span.text == 'hg19':
                        ref_hg19 = span.attrs['title']
                    elif span.text == 'hg38':
                        ref_hg38 = span.attrs['title']

                rsid = '.'
                if tds[-1].select('a[target="dbSNP"]'):
                    rsid = tds[-1].select('a[target="dbSNP"]')[0].attrs['href'].split('rs=')[-1]

                context = {}
                for each in 'hgmd_id mutation phenotype display reports pmids ref_hg19 ref_hg38 rsid'.split():
                    context.update({
                        each: locals()[each]
                    })
                contexts.append(context)

        return contexts

    @try_again()
    def get_reports(self, hgmd_id, session, proxies):

        # print '\033[32m> get_report: {} \033[0m'.format(hgmd_id)

        soup = self.get_soup(MUTATION_URL + hgmd_id, session, proxies)
        report_list = []
        pmid_list = []
        for tr in soup.select('table.gene')[2].select('tr.odd'):

            rc1 = re.compile(
                r'(?:^\d+. )|(?: PubMed:)'
            )
            rc2 = re.compile(
                r'(?:^\d+.\s)|(?:no PubMed ID$)|(?: -  LSDB$)|(?: Springer$)|(?:\(meeting abstract\)$)'
            )
            rc2 = re.compile(
                r'(?:^\d+.\s)|(?:no PubMed ID$)|(?: -  LSDB$)|(?: Springer$)|(?:\(meeting abstract\)$)'
            )
            content = tr.select('td')[0].text.strip()
            if 'PubMed:' in content:
                result = rc1.split(content)
            else:
                result = rc2.split(content)

            if len(result) != 3:
                print repr(content)
                print result
                raise Exception('parse pmid failed for ' + hgmd_id)

            cit = result[1].strip()
            report_list.append(cit)

            pmid = result[2].strip()
            if pmid in ('ID', ''):
                pmid = '.'
            pmid_list.append(pmid)
        return report_list, pmid_list

    @try_again(10)
    def get_soup(self, url, session, proxies, method='GET', **kwargs):
        # 每次请求都重新登录
        # session, proxies = self.login()
        # print '\033[32m> get_soup: {} \033[0m'.format(url)

        timeout = 15
        if method == 'GET':
            response = session.get(
                url, timeout=timeout, proxies=proxies, **kwargs)
        elif method == 'POST':
            response = session.post(
                url, timeout=timeout, proxies=proxies, **kwargs)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        return soup

    @try_again()
    def login(self):

        proxies = self.proxy.get_random_proxies(protocol='https')

        data = {'password': PASSWORD, 'email': EMAIL}
        session = requests.session()
        timeout = 10
        response = session.post(LOGIN_URL, data=data, headers=HEADERS, proxies=proxies, timeout=timeout)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        if 'Login successful' in soup.select('table td[align="center"]')[0].text:
            print '\033[32m> login successful with proxies: {}\033[0m'.format(proxies)
            return session, proxies
        else:
            print '\033[31m> login failed with proxies: {}\033[0m'.format(proxies)

    @property
    def genelist(self):
        if args['gene']:
            genelist = self.args['gene'].split(',')
        else:
            if os.path.isfile(self.args['genelist']):
                print 'use local genelist file: {}'.format(
                    self.args['genelist'])
                with open(self.args['genelist']) as f:
                    genelist = [each for each in f.read().strip().split('\n')]
            elif self.args['craw_genelist']:
                print 'crawling all genelist from HGMD...'
                genelist = self.craw_genelist(self.args['genelist'])
            else:
                print 'no genelist was supplied, you can craw all genelist from hgmd with -cg argument'
                exit(0)
        print 'the length of genelist: {}'.format(len(genelist))
        return genelist

    @try_again()
    def craw_genelist(self, outfile):

        limit = '6000'

        genelist = []
        with open(outfile, 'w') as out:
            for display in ('DM', 'P', 'DMP'):
                print 'crawling {} gene...'.format(display)
                data = {'display': display, 'limit': limit}

                session, proxies = self.login()
                soup = self.get_soup(GENE_TAG_URL, session, proxies, method='POST', data=data)

                for tr in soup.select('table')[1].select('tr')[1:]:
                    tds = tr.select('td')
                    gene_symbol = tds[0].select('a')[0].text
                    genelist.append(gene_symbol)
                    out.write(gene_symbol + '\n')
                    print 'find gene: {}'.format(gene_symbol)

        return genelist


def main():

    proxy = Proxy(ipPage=30)
    hgmd = CrawHGMD(args, proxy)
    hgmd.start()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='craw all gene information from HGMD',
        prog='craw_gene',
        epilog='contact: suqingdong@novogene.com',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gene', help='the input gene')
    parser.add_argument(
        '-gl',
        '--genelist',
        help=
        'the input genelist file, one gene per line[default="%(default)s"]',
        default='HGMD.genelist')
    parser.add_argument(
        '-o',
        '--out',
        help='the prefix of outfile[default=%(default)s]',
        default='HGMD')
    parser.add_argument(
        '-cg',
        '--craw-genelist',
        action='store_true',
        help='Craw all genelist from HGMD')

    args = vars(parser.parse_args())

    main()
