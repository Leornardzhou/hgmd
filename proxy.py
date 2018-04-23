#!/usr/bin/env python
#-*- enconding: utf-8 -*-
'''
Usage:
    from proxy import Proxy
    proxy = Proxy(maxTry=10, ipPage=10, linkSpeed=1)
    soup = proxy.get_soup(your_url, protocol='http')
'''

import os
import time
import random
import bs4
import requests

import json
from collections import defaultdict


class Proxy:

    def __init__(self, maxTry=5, linkSpeed=2, ipPage=1, IP_POOL_JSON='ip_pool.json'):
        self.headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36'}
        self.maxTry = maxTry
        self.linkSpeed = linkSpeed
        self.ipPage = ipPage
        self.IP_POOL_JSON = IP_POOL_JSON

        self.ip_pool = self.generate_ip_pool()  # must be the last attr, only generate when initial an instance

    def generate_ip_pool(self):
        '''
        Generate a proxy IP Pool when initial.
        '''
        if os.path.exists(self.IP_POOL_JSON):
            with open(self.IP_POOL_JSON) as f:
                ip_pool = json.load(f)
                print( 'Use local ip pool: ' + self.IP_POOL_JSON )
                print( 'HTTP IP Pools Counts: ', len(ip_pool.get('http')) )
                print( 'HTTPS IP Pools Counts: ', len(ip_pool.get('https')) )
                return ip_pool

        ip_pool = defaultdict(list)
        for page in range(1, self.ipPage+1):
            ip_url = 'http://www.xicidaili.com/nn/' + str(page)
            response = requests.get(ip_url, headers=self.headers)
            soup = bs4.BeautifulSoup(response.text, 'html.parser')
            ips = soup.select('#ip_list tr')            
            for ip_info in ips[1:]:
                tds = ip_info.select('td')
                ip = tds[1].text
                port = tds[2].text
                protocol = tds[5].text.lower()
                speed = float(tds[6].findChild().attrs['title'][:3])
                if speed < self.linkSpeed:
                    proxies = { protocol: '{}://{}:{}'.format(protocol, ip, port) }
                    ip_pool[protocol].append(proxies)
        print( 'HTTP IP Pools Counts: ', len(ip_pool.get('http')) )
        print( 'HTTPS IP Pools Counts: ', len(ip_pool.get('https')) )

        with open(self.IP_POOL_JSON, 'w') as f:
            json.dump(ip_pool, f)

        return ip_pool

    def get_random_proxies(self, protocol='http'):
        '''
        Get a random proxies for requests.get.
        '''
        proxy_list = self.ip_pool.get(protocol)
        random_proxies = random.choice(proxy_list)
        return random_proxies

    def get_soup(self, url, enconding='utf-8', protocol='http'):
        '''
        Requst with a random proxies, and return a soup.
        '''
        for n in range(1, self.maxTry+1):
            try:
                proxies = self.get_random_proxies(protocol=protocol)
                print('> Try connect with proxies: ', proxies)
                response = requests.get(url, headers=self.headers, proxies=proxies)
                soup = bs4.BeautifulSoup(response.content.decode(enconding), 'html.parser')
                print('>>> Connection successful with IP: ' + proxies.values()[0])
                return soup
            except Exception as e:
                print('* Connection failed {}st time: {}'.format(n, e))
        print('Opps! Fail too many times(%s)! Yon can retry later or change the number of maxTry.' % self.maxTry)


if __name__ == '__main__':

    proxy = Proxy(maxTry=5, ipPage=40)

    # Test NCBI
    url = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=rs371061770'
    soup = proxy.get_soup(url)
    print 'Result: ' + soup.select('.clinsig_confidence')[0].text

    # Test GeneCards
    url = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=BTK'
    soup = proxy.get_soup(url, protocol='https')
    print 'Result: ' + soup.select('#summaries')[0].text.strip()
