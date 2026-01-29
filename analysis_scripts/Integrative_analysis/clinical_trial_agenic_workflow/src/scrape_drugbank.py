#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import undetected_chromedriver as uc
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException, ElementClickInterceptedException
from bs4 import BeautifulSoup
import pandas as pd
import time


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--output',type=str,dest="output",help="output file",default=None,required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def main():
    global args
    args = parse_arg()

    # instantiate a Chrome browser
    url = 'https://go.drugbank.com/conditions/DBCOND0049114/'  # Your URL
    driver = uc.Chrome(
        use_subprocess=False,
        headless=True,
    )

    # visit the target page
    driver.get(url)

    # wait for the interstitial page to load
    time.sleep(10)

    # take a screenshot of the current page and save it
    driver.save_screenshot("cloudflare-challenge.png")
    exit(1)

    # close the browser
    driver.quit()
    exit(1)

    options = webdriver.ChromeOptions()
    options.add_argument("--headless")   # example
    options.add_argument("--user-data-dir=/home/yangz6/.config/google-chrome/Default") # example
    options.add_argument("--disable-blink-features=AutomationControlled")
    driver = webdriver.Chrome(options=options)
    url = 'https://go.drugbank.com/conditions/DBCOND0049114/'  # Your URL
    driver.get(url)
    # check if the page is loaded
    time.sleep(12)  # Wait for the page to load
    print(driver.title)
    exit(1)
    # 
    wait = WebDriverWait(driver, 10)
    table_element = wait.until(EC.presence_of_element_located((By.ID, 'DataTables_Table_1')))
    time.sleep(2)

    # Initial parse
    html = driver.page_source
    soup = BeautifulSoup(html, 'html.parser')
    table = soup.find('table', id='DataTables_Table_1')
    # 
    all_rows = []
    for row in table.find_all('tr'):
        cols = [col.text.strip() for col in row.find_all(['td', 'th'])]
        if cols:
            all_rows.append(cols)
    print(all_rows)




    
if __name__=="__main__":
    main()

