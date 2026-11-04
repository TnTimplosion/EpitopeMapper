#!/usr/bin/env python3
"""Serve the app over HTTP for local use.

    python tools/serve.py            # http://localhost:8765
    python tools/serve.py -p 9000

The app itself runs straight from the file system, so this is only needed for the pages in tests/,
which read their fixtures with fetch() and therefore need a real origin. Responses carry
Cache-Control: no-store so an edited file shows up on the next reload rather than being served from
the browser's memory cache.
"""
import argparse
import http.server
import pathlib
import socketserver

ROOT = pathlib.Path(__file__).resolve().parent.parent


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, must-revalidate')
        super().end_headers()

    def log_message(self, fmt, *args):
        if '200' not in (args[1] if len(args) > 1 else ''):
            super().log_message(fmt, *args)


class Server(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-p', '--port', type=int, default=8765)
    args = ap.parse_args()
    with Server(('127.0.0.1', args.port), Handler) as httpd:
        print('serving %s at http://localhost:%d' % (ROOT, args.port))
        print('  app    http://localhost:%d/index.html' % args.port)
        print('  tests  http://localhost:%d/tests/core_test.html' % args.port)
        print('         http://localhost:%d/tests/accuracy_test.html' % args.port)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print('\nstopped')


if __name__ == '__main__':
    main()
