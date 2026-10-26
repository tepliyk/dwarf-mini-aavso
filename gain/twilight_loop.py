# -*- coding: utf-8 -*-
"""Сутінкова серія для STWG (замір gain): кожні INTERVAL с запускає RAW-серію через dwarf_ops.capture_raw
з фіксованою витримкою/gain, до часу STOP_AT. Лог у twilight_loop.log. Запуск: python -X utf8 twilight_loop.py"""
import asyncio, datetime as dt, sys, os
sys.path.insert(0, r'C:\astro\dwarf_mcp')
import dwarf_ops as ops

COUNT, EXP, GAIN, MODE, IR = 60, '1/4', 0, 10, 1
INTERVAL = 180          # с між стартами серій
STOP_AT = dt.time(17, 22)
LOG = r'C:\astro\stwg\twilight_loop.log'

def log(s):
    line = '%s %s' % (dt.datetime.now().strftime('%H:%M:%S'), s)
    print(line); open(LOG, 'a', encoding='utf-8').write(line + '\n')

async def main():
    n = 0
    while dt.datetime.now().time() < STOP_AT:
        t0 = dt.datetime.now()
        try:
            r = await ops.capture_raw(COUNT, EXP, GAIN, MODE, IR)
            n += 1
            log('серія %d: ok=%s folder=%s files=%s' % (n, r.get('ok'), r.get('folder'), r.get('files')))
        except Exception as e:
            log('ПОМИЛКА: %r' % e)
            try:
                await ops.stop()
            except Exception:
                pass
        wait = INTERVAL - (dt.datetime.now() - t0).total_seconds()
        if wait > 0:
            await asyncio.sleep(wait)
    log('стоп за часом, серій %d' % n)

asyncio.run(main())
