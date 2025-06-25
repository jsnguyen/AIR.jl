import marimo

__generated_with = "0.12.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    from pathlib import Path
    import shutil

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.table import Table, Column
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    from pykoa.koa import Koa
    return Column, Koa, Path, SkyCoord, Table, fits, mo, np, plt, shutil, u


@app.cell
def _(__file__):
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)
    cwd = os.getcwd()
    cwd
    return cwd, dir_path, os


@app.cell
def _():
    dates_filename = 'AS_209_KOA_dates.txt'

    dates = []
    with open(dates_filename, 'r') as f:
        lines = f.readlines()
        for d in lines:
            dates.append(d.strip())

    dates
    return d, dates, dates_filename, f, lines


@app.cell
def _(SkyCoord):
    coord = SkyCoord.from_name("AS 209", parse=True)
    return (coord,)


@app.cell
def _(Koa, coord):
    table_filename = 'nirc2_as_209.vot'
    pos_query = f'circle {coord.ra.value} {coord.dec.value} 0.5'
    Koa.query_position('nirc2', pos_query, table_filename, overwrite=True, format='votable')
    return pos_query, table_filename


@app.cell
def _(Table, table_filename):
    rec = Table.read(table_filename,format='votable')
    koa_dates = sorted(set(rec['date_obs'].data))
    print(koa_dates)
    #print(dates)
    return koa_dates, rec


@app.cell
def _(Koa, koa_dates):
    for date in koa_dates:
        date_table_filename = f'{date}_nirc2_as_209.vot'
        Koa.query_date('nirc2', date, date_table_filename, overwrite=True, format='votable' )
    return date, date_table_filename


@app.cell
def _(Table, koa_dates):
    def _():
        for date in koa_dates:
            date_table_filename = f'{date}_nirc2_as_209.vot'
            rec = Table.read(date_table_filename, format='votable')
            print(rec.length)

    _()
    return


@app.cell
def _(Koa, koa_dates):
    def _():
        for date in koa_dates:
            date_table_filename = f'{date}_nirc2_as_209.vot'
            Koa.download(date_table_filename, 'votable', 'nirc2_data/'+date, calibfile=1, calibdir=0)

    _()
    return


@app.cell
def _(Path, fits, shutil):
    def _():
        data_folder = Path('nirc2_data')

        for date in sorted(data_folder.glob('*/')):
            for filepath in sorted((date / 'lev0').glob('*.fits')):
                header = fits.getheader(filepath)
                #print(header['FILENAME'])
                new_folder = filepath.parents[1]/'raw'
                new_folder.mkdir(parents=True, exist_ok=True)
                new_filepath = new_folder / header['FILENAME']
                if new_filepath.exists():
                    print('EXISTS, RENAMING')
                    new_filepath = new_folder / header['FILENAME'].replace('.fits', '_dup.fits')
                    #continue
                print(f'Moving {filepath} -> {new_filepath}')
                shutil.move(filepath, new_filepath)
            
        return


    _()
    return


@app.cell
def _(Path, shutil):
    def _():
        data_folder = Path('nirc2_data')

        for folder in sorted(data_folder.glob('*/lev0/')):
            if list(folder.glob('.*fits')) == []:
                print(f'REMOVING {folder}')
                shutil.rmtree(folder)
            
        return


    _()
    return


if __name__ == "__main__":
    app.run()
