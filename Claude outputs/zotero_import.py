#!/usr/bin/env python3
"""Create the Pd-membrane WGS collection tree in Zotero and add the items.

Standard library only. Default is a DRY RUN that prints the plan and makes no
network call to your library. To write to Zotero:

    python zotero_import.py --user-id 1234567 --api-key YOUR_KEY --apply

Get the numeric user ID and a key with write access at
https://www.zotero.org/settings/keys  (Create new private key -> allow write access).
Run it once. A second run stops if the root collection already exists.
"""
import argparse, json, os, sys, time, uuid, urllib.request, urllib.error
API = 'https://api.zotero.org'  # do not edit; the key goes on the command line or in ZOTERO_API_KEY

def call(method, path, key, body=None, extra=None, retries=4):
    h = {'Zotero-API-Version': '3', 'Content-Type': 'application/json'}
    if key: h['Zotero-API-Key'] = key
    if extra: h.update(extra)
    data = json.dumps(body).encode() if body is not None else None
    for k in range(retries):
        try:
            req = urllib.request.Request(API + path, data=data, headers=h, method=method)
            with urllib.request.urlopen(req, timeout=60) as r:
                bo = r.headers.get('Backoff')
                if bo: time.sleep(float(bo))
                txt = r.read().decode()
                return json.loads(txt) if txt else {}
        except ValueError as e:
            sys.exit(f'Bad request URL ({e}). Check that the script file is unedited.')
        except urllib.error.HTTPError as e:
            if e.code in (429, 503):
                time.sleep(float(e.headers.get('Retry-After', 5)) + 1); continue
            sys.exit(f'HTTP {e.code} on {method} {path}: {e.read().decode()[:400]}')
        except urllib.error.URLError as e:
            sys.exit(f'Network error on {method} {path}: {e.reason}')
    sys.exit('Gave up after repeated rate-limit responses.')

def templates_ok(items):
    """Compare field names with Zotero's own item templates (no key needed)."""
    bad = {}
    cache = {}
    for it in items:
        t = it['item']['itemType']
        if t not in cache:
            cache[t] = set(call('GET', f'/items/new?itemType={t}', None).keys())
        extra = [k for k in it['item'] if k not in cache[t] and k not in ('collections',)]
        if extra: bad.setdefault((t, tuple(sorted(extra))), []).append(it['review_id'])
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--user-id'); ap.add_argument('--api-key')
    ap.add_argument('--data', default='corpus_items.json')
    ap.add_argument('--apply', action='store_true', help='write to your library')
    ap.add_argument('--skip-template-check', action='store_true')
    a = ap.parse_args()
    a.api_key = a.api_key or os.environ.get('ZOTERO_API_KEY')
    if a.apply and a.api_key and a.api_key.upper().startswith('YOUR_'):
        sys.exit('The key is still the placeholder text. Use the real key.')
    d = json.load(open(a.data, encoding='utf-8'))
    items, paths = d['items'], d['collections']
    print(f"Plan: 1 root + {len(paths)} sub-collections, {len(items)} items.")
    for p in paths: print('  ', p)
    if not a.apply:
        print('\nDry run only. Add --user-id, --api-key and --apply to create them.'); return
    if not (a.user_id and a.api_key): sys.exit('--user-id and --api-key are required with --apply')
    if not a.skip_template_check:
        bad = templates_ok(items)
        if bad:
            print('Field names not in the Zotero template (fix or re-run with --skip-template-check):')
            for k, v in bad.items(): print('  ', k, v[:5])
            sys.exit(1)
        print('Template check passed.')
    base = f'/users/{a.user_id}'
    root = d['root']
    top = call('GET', base + '/collections/top?limit=100', a.api_key)
    if any(c['data']['name'] == root for c in top): sys.exit(f'Root collection "{root}" already exists. Stopping.')
    wt = lambda: {'Zotero-Write-Token': uuid.uuid4().hex}
    r = call('POST', base + '/collections', a.api_key, [{'name': root, 'parentCollection': False}], wt())
    rootkey = r['success']['0']; print('Root collection key', rootkey)
    subs = [p.split('/', 1)[1] for p in paths]
    r = call('POST', base + '/collections', a.api_key, [{'name': s, 'parentCollection': rootkey} for s in subs], wt())
    keymap = {f'{root}/{subs[int(i)]}': k for i, k in r['success'].items()}
    ok = fail = 0
    for i in range(0, len(items), 50):
        chunk = items[i:i + 50]; payload = []
        for it in chunk:
            o = dict(it['item']); o['collections'] = [keymap[c] for c in it['collections']]
            payload.append(o)
        r = call('POST', base + '/items', a.api_key, payload, wt())
        ok += len(r.get('success', {})); fail += len(r.get('failed', {}))
        for idx, f in r.get('failed', {}).items(): print('FAILED', chunk[int(idx)]['review_id'], f)
    print(f'Done: {ok} items created, {fail} failed.')

if __name__ == '__main__': main()
