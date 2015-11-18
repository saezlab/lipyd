
#!/bin/bash

import pycurl
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
try:
    import cPickle as pickle
except:
    import pickle
import sys
import os
import re
import urllib
import urlparse
import gzip
import zipfile
import tarfile
import hashlib

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

show_cache = False

def is_quoted(string):
    '''
    From http://stackoverflow.com/questions/1637762/test-if-string-is-url-encoded-in-php
    '''
    test = string
    while(urllib.unquote(test) != test):
        test = urllib.unquote(test)
    return urllib.quote(test, '/%') == string or urllib.quote(test) == string

def is_quoted_plus(string):
    test = string
    while(urllib.unquote_plus(test) != test):
        test = urllib.unquote_plus(test)
    return urllib.quote_plus(test, '&=') == string or urllib.quote_plus(test) == string

def url_fix(s, charset='utf-8', force = False):
    """
    From http://stackoverflow.com/a/121017/854988
    """
    if isinstance(s, unicode):
        s = s.encode(charset, 'ignore')
    scheme, netloc, path, qs, anchor = urlparse.urlsplit(s)
    if force or not is_quoted(path):
        path = urllib.quote(path, '/%')
    if force or not is_quoted_plus(qs):
        qs = urllib.quote_plus(qs, '&=')
    return urlparse.urlunsplit((scheme, netloc, path, qs, anchor))

def print_debug_info(debug_type, debug_msg, truncate = 1000):
    sys.stdout.write("debug(%d): %s\n" % (debug_type, debug_msg[:truncate]))
    sys.stdout.flush()

#class Dataio(object):
#    
#    __init__(self,mapper=None):

def get_headers(header_list):
    headers = {}
    for header_line in header_list:
        if ':' not in header_line:
            continue
        name, value = header_line.split(':', 1)
        name = name.strip()
        value = value.strip()
        name = name.lower()
        headers[name] = value
    return headers

def get_jsessionid(headers):
    rejsess = re.compile(r'.*(JSESSIONID=[A-Z0-9]*)')
    for hdr in headers:
        jsess = rejsess.findall(hdr)
        if len(jsess) > 0:
            return ['Cookie: %s'%jsess[0]]

def get_xsessionid(headers):
    pass

def curl(url, silent = True, post = None, req_headers = None, cache = True, 
        debug = False, outf = None, compr = None, encoding = None, 
        files_needed = None, timeout = 300, init_url = None, 
        init_fun = 'get_jsessionid', follow = True, large = False,
        override_post = False, init_headers = False, 
        write_cache = True, force_quote = False):
    url = url_fix(url, force = force_quote)
    if init_url is not None:
        init_url = url_fix(init_url, force = force_quote)
    # either from cache or from download, we load the data into StringIO:
    multifile = False
    domain = url.replace('https://', '').replace('http://','').\
        replace('ftp://','').split('/')[0]
    # first try to find file in cache:
    if cache or write_cache:
        # outf param is to give a unique name to data
        # downloaded previously by post requests
        outf = outf if outf is not None else url.split('/')[-1].split('?')[0]
        poststr = '' if post is None else \
            '?' + '&'.join(sorted([i[0]+'='+i[1] for i in post.items()]))
        try:
            urlmd5 = hashlib.md5(url+poststr).hexdigest()
        except UnicodeEncodeError:
            urlmd5 = hashlib.md5(('%s%s' % (url, poststr)).encode('utf-8')).hexdigest()
        if not os.path.exists(os.path.join(os.getcwd(),'cache')):
            os.mkdir(os.path.join(os.getcwd(),'cache'))
        cachefile = os.path.join(os.getcwd(),'cache',urlmd5+'-'+outf)
        if show_cache:
            sys.stdout.write('\tFor URL %s\n' % url)
            sys.stdout.write('\tChache file is %s' % cachefile)
        usecache = True if os.path.exists(cachefile) and cache else False
        # load from cache:
        if usecache:
            if not silent:
                sys.stdout.write('\t:: Loading %s from cache, previously '\
                    'downloaded from %s\n'%(outf,domain))
                sys.stdout.flush()
            if large:
                result = open(cachefile, 'rb')
            else:
                with open(cachefile,'rb') as f:
                    result = StringIO()
                    result.write(f.read())
    else:
        usecache = False
    # if not found in cache, download with curl:
    if not usecache:
        headers = []
        if not init_url and large:
            result = open(cachefile, 'wb')
        else:
            result = StringIO()
        c = pycurl.Curl()
        if init_url:
            c.setopt(c.URL, init_url)
        else:
            try:
                c.setopt(c.URL, url)
            except:
                return url
        c.setopt(c.FOLLOWLOCATION, follow)
        c.setopt(c.CONNECTTIMEOUT, 15)
        c.setopt(c.TIMEOUT, timeout)
        if override_post:
            if req_headers is None: req_headers = []
            req_headers.append('X-HTTP-Method-Override: GET')
        if type(req_headers) is list:
            c.setopt(c.HTTPHEADER, req_headers)
        c.setopt(c.WRITEFUNCTION, result.write)
        c.setopt(c.HEADERFUNCTION, headers.append)
        # if debug is necessary:
        if debug:
            c.setopt(pycurl.VERBOSE, 1)
            c.setopt(pycurl.DEBUGFUNCTION, print_debug_info)
        if type(post) is dict:
            postfields = urllib.urlencode(post)
            c.setopt(c.POSTFIELDS, postfields)
            c.setopt(c.POST, 1)
        if not silent:
            sys.stdout.write('\t:: Downloading data from %s. Waiting for reply...' % \
                domain)
            sys.stdout.flush()
        for i in xrange(3):
            try:
                if debug:
                    sys.stdout.write('\t:: bioigraph.dataio.curl() :: attempt #%u\n' % i)
                    sys.stdout.flush()
                c.perform()
                if url.startswith('http'):
                    status = c.getinfo(pycurl.HTTP_CODE)
                    if status == 200:
                        break
                if url.startswith('ftp'):
                    status = 500
                    for h in headers:
                        if h.startswith('226'):
                            status = 200
                            break
            except pycurl.error as (errno, strerror):
                status = 500
                sys.stdout.write('\tPycURL error: %u, %s\n' % (errno, strerror))
                sys.stdout.flush()
        c.close()
    # sometimes authentication or cookies are needed to access the target url:
    if init_url and not usecache:
        if not silent:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + 'Success.\n')
            sys.stdout.flush()
        # here, you may define a custom function to fetch 
        # the authentication data from cookies/headers, 
        # and return with headers for the main request:
        req_headers = globals()[init_fun](headers)
        if init_headers: return req_headers
        return curl(url = url, req_headers = req_headers, silent = silent, 
            debug = debug, outf = outf, compr = compr, encoding = encoding, 
            files_needed = files_needed, timeout = timeout, large = large,
            write_cache = write_cache)
    # get the data from the file downloaded/loaded from cache:
    if usecache or status == 200:
        if type(result) is file:
            fname = result.name
            result.close()
            result = open(fname, 'r')
        # find out the encoding:
        if encoding is None:
            if not usecache:
                headers = get_headers(headers)
                encoding = None
                if 'content-type' in headers:
                    content_type = headers['content-type'].lower()
                    match = re.search('charset=(\S+)', content_type)
                    if match:
                        encoding = match.group(1)
                if encoding is None:
                    if url.startswith('ftp'):
                        encoding = 'utf-8'
                    else:
                        encoding = 'iso-8859-1'
            else:
                # in case of using the cache:
                encoding = 'utf-8'
        if not silent and not usecache:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + 'Success.\n')
            sys.stdout.flush()
        result.seek(0)
        if url.endswith('tar.gz') or url.endswith('tgz') or compr == 'tgz':
            multifile = True
            results = {}
            res = tarfile.open(fileobj = result, mode = 'r:gz')
            membs = res.getmembers()
            for m in membs:
                if (files_needed is None or m.name in files_needed) \
                    and m.size != 0:
                    # m.size is 0 for dierctories
                    this_file = res.extractfile(m)
                    if large:
                        results[m.name] = this_file
                    else:
                        results[m.name] = this_file.read()
                        this_file.close()
            if not large:
                res.close()
        elif url.endswith('gz') or compr == 'gz':
            res = gzip.GzipFile(fileobj=result, mode='rb')
            if not large:
                res = res.read()
                try:
                    res = res.decode(encoding)
                    res = res.encode('utf-8')
                except:
                    # better to proceed even if there is some trouble with encodings...
                    pass
        elif url.endswith('zip') or compr == 'zip':
            multifile = True
            results = {}
            res = zipfile.ZipFile(result,'r')
            membs = res.namelist()
            for m in membs:
                if files_needed is None or m in files_needed:
                    this_file = res.open(m)
                    if large:
                        results[m] = this_file
                    else:
                        results[m] = this_file.read()
                        this_file.close()
            res.close()
        else:
            if large:
                res = result
            else:
                res = result.getvalue()
        if not multifile:
            results = {'one': res}
        if not large:
            for k in results.keys():
                # handle files with CR line endings:
                if '\r' in results[k] and '\n' not in results[k]:
                    results[k] = results[k].replace('\r','\n')
                else:
                    results[k] = results[k].replace('\r','')
                if 'encoding' != 'utf-8':
                    try:
                        results[k] = results[k].decode(encoding).encode('utf-8')
                    except:
                        pass
        if (cache or write_cache) and not usecache and not large:
            for k in results.keys():
                if not multifile and not url.endswith('gz') and compr != 'gz':
                # write the decoded data back to StringIO
                    result.truncate(0)
                    result.write(results[k])
                # if cache is turned on, but data is not from cache,
                # place it there to make available next time:
                result.seek(0)
                with open(cachefile,'wb') as f:
                    f.write(result.getvalue())
        res = results if multifile else results['one']
    else:
        # download error:
        if not silent:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + \
                'Failed. (Status: %u)\n'%status)
            if status > 200:
                sys.stdout.write('\t# URL: %s\n\t# POST: %s\n' % \
                    (url, '' if type(post) is not dict else urllib.urlencode(post)))
            sys.stdout.flush()
        res = None
    # returns raw data, dict of file names and raw data in case of 
    # multiple file archives, or file object in case of large files:
    return res
