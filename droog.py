"""Generate drug name misspellings"""
#
# Copyright 2017 Robert D. Hogan
# rhogan@terminologix.com
#
# This file is part of Droog
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import  sys, os, argparse
import  matplotlib.pyplot as plt
from    noisy import *

def main():
    """CLI for generating drug spelling variants"""
    parser = argparse.ArgumentParser(description='Drug misspelling generator')
    parser.add_argument('-d', '--drug', help='correctly spelled drug name', \
                        required=True)
    parser.add_argument('-g', '--generate', action='store_true', \
                        help='generate variants with probabilities')
    parser.add_argument('-f', '--filter', action='store_true', \
                        help='filter variants with Google CSE')
    parser.add_argument('-a', '--analyze', action='store_true', \
                        help='print table of terms vs. cumulative counts')
    parser.add_argument('-l', '--list', type=list_arg_type, default=False, \
                        help='list n, n%% or \'all\' filtered terms')
    parser.add_argument('-p', '--plot', action='store_true', \
                        help='plot cumulative page count vs. list number')

    # parse command line parameters
    args = vars(parser.parse_args())
    drug = args['drug']
    path = '.\\results\\' + drug

    # error check parameters
    optionals = ['generate', 'filter', 'plot', 'analyze', 'list']
    existing_drugs = os.listdir('.\\results')
    arg_count = sum([1 for o in optionals if args[o] != False])
    if arg_count != 1:
        error_exit('must have one and only one optional argument')
    if not args['generate'] and drug not in existing_drugs:
        error_exit('must generate variants first')
    for option in optionals:
        if args[option] != False:
            command = option

    # execute 'generate' command
    if command == 'generate':
        if drug not in existing_drugs:
            os.mkdir(path)
        else:
            verify_files(path, drug, 'd1', 'd2', 'likely', \
                         confirm_overwrite=True)
        d1, d2 = generate(drug)
        save_candidates(path, drug, 'd1', d1)
        save_candidates(path, drug, 'd2', d2)
        likely = d1[:int(len(d1)*D1_LIKELY_PCT/100)] + \
                 d2[:int(len(d2)*D2_LIKELY_PCT/100)]
        save_candidates(path, drug, 'likely', likely, term_only=True)

    # execute 'filter' command
    if command == 'filter':
        verify_files(path, drug, 'filtered', confirm_overwrite=True)
        verify_files(path, drug, 'd1', 'd2', confirm_exist=True)
        d1, d2 = read_candidates(path, drug)
        filtered = filter_variants(drug, d1, d2, limit=0, progress=True)
        with open(make_file_name(path, drug, 'filtered'), 'w') as file:
            for row in filtered:
                file.write(row[0] + ',' + str(row[1]) + '\n')

    # execute 'analyze' command
    if command == 'analyze':
        correct_pages, total_pages, total_misspells, table = \
            analyze_pages(path, drug)
        # print summary
        print('\n{:<28} {:<15}'.format('drug:', drug))
        print('{:<28} {:<5d}'.format('total # of pages:', total_pages))
        print('{:<28} {:d} ({:.1f}%)'.format('pages with correct spelling:', \
                                            correct_pages, \
                                            correct_pages/total_pages*100))
        print('{:<28} {:d} ({:.1f}%)'.format('pages with misspelling(s):', \
                                            total_misspells, \
                                            total_misspells/total_pages*100))
        print('{:<28} {:<15}'.format('# of unique misspellings:', len(table)))
        print()
        # print table header
        print('{}{}{}{}'.format(' '*32, 'all pages', ' '*11, \
                                'misspelled pages'))
        print(' #  misspelling     pages      %       cum %         %       ' \
              '    cum %')
        print('--- -----------     -----    ----      -----       -----     ' \
              '    -----')
        # print table
        template = '{:>3d} {:<15} {:>4d}   {:>5.1f}%     {:>5.1f}%      ' \
                   '{:>5.1f}%        {:>5.1f}%'
        for i, row in enumerate(table):
            print(template.format(i+1, row['misspelling'], row['pages'], \
                                  row['pct_pages'], row['cum_pct_pages'], \
                                  row['pct_misspells'], row['cum_pct_misspells']))

    # execute 'list' command
    if command == 'list':
        verify_files(path, drug, 'filtered', confirm_exist=True)
        filtered = read_delimited_list(path, drug + '_filtered_variants.txt')
        total_misspells = sum([int(x[1]) for x in filtered if x[0] != drug])
        n_type, n = args[command]
        misspell_pages = 0
        for i, row in enumerate(filtered):
            if row[0] != drug:
                misspell_pages += int(row[1])
                if n_type == 'all' or \
                           n_type == 'num' and \
                           i <= n or n_type == '%' and \
                           misspell_pages/total_misspells <= n/100:
                    print(row[0])

    # execute 'plot' command
    if command == 'plot':
        correct_pages, total_pages, total_misspells, table = \
            analyze_pages(path, drug)
        fig, ax = plt.subplots(1, 2)
        plt.subplots_adjust(wspace=0.3)
        fig.suptitle('Misspelling Capture Curves', fontsize=18)
        x = [i for i in range(len(table))]
        y1 = [row['cum_pct_pages'] for row in table]
        y2 = [row['cum_pct_misspells'] for row in table]
        ax[0].set_autoscaley_on(False)
        ax[0].set_ylim([0, 100])
        ax[0].plot(x, y1, color='r', lw=2)
        ax[0].set_xlabel('# of Misspellings')
        ax[0].set_ylabel('Cumulative % of All Pages')
        ax[0].text(0.1, 0.1, 'includes correct spelling', ha='left', \
                   va='center', transform=ax[0].transAxes)
        ax[1].set_autoscaley_on(False)
        ax[1].set_ylim([0, 100])
        ax[1].plot(x, y2, color='b', lw=2)
        ax[1].set_xlabel('# of Misspellings')
        ax[1].set_ylabel('Cumulative % of Misspelled Pages')
        plt.show()

def analyze_pages(path, drug):
    """create statistics table for drug name misspellings"""
    verify_files(path, drug, 'filtered', confirm_exist=True)
    filtered = read_delimited_list(path, drug + '_filtered_variants.txt')
    filtered = [x for x in filtered if x[0][0] != '#']
    # compute totals
    correct_pages = int(filtered[0][1])
    total_pages = sum([int(x[1]) for x in filtered])
    total_misspells = sum([int(x[1]) for x in filtered if x[0] != drug])
    # generate table
    sum_pages, sum_misspells, table = correct_pages, 0, []
    for row in filtered:
        misspelling, pages, new_row = row[0], int(row[1]), {}
        if misspelling != drug:
            sum_pages += pages
            sum_misspells += pages
            new_row['misspelling'] = misspelling
            new_row['pages'] = pages
            new_row['pct_pages'] = pages/total_pages*100
            new_row['cum_pct_pages'] = sum_pages/total_pages*100
            new_row['pct_misspells'] = pages/total_misspells*100
            new_row['cum_pct_misspells'] = sum_misspells/total_misspells*100
            table.append(new_row)
    return correct_pages, total_pages, total_misspells, table

def verify_files(path, drug, *args, **kwargs):
    """check if any files exist and get permission to overwrite"""
    prefix = path + '\\' + drug + '_'
    postfix = '_variants.txt'
    files = [prefix + x + postfix for x in args]
    count = sum(1 for x in files if os.path.isfile(x))
    if 'confirm_overwrite' in kwargs:
        if count > 0:
            answer = input("overwrite files? (y/n): ")
            if answer != 'y':
                exit()
    if 'confirm_exist' in kwargs:
        if count != len(files):
            error_exit('required files do not exist')
    return

def save_candidates(path, drug, v_type, candidates, term_only=False):
    """save spelling variants as comma-delimited list"""
    file_name = make_file_name(path, drug, v_type)
    with open(file_name, 'w') as file:
        for row in candidates:
            if term_only:
                file.write(row[0] + '\n')
            else:
                file.write(','.join(str(e) for e in row) + '\n')

def read_candidates(path, drug):
    """read d1 and d2 candidate files"""
    d1_file = make_file_name(path, drug, 'd1')
    d2_file = make_file_name(path, drug, 'd2')
    d1, d2 = [], []
    with open(d1_file, 'r') as file:
        for line in file:
            d1.append(line.split(','))
    with open(d2_file, 'r') as file:
        for line in file:
            d2.append(line.split(','))
    return d1, d2

def make_file_name(path, drug, file_type):
    """generate full file name with path for a drug"""
    return path + '\\' + drug + '_' + file_type + '_variants.txt'

def list_arg_type(string):
    """type check argument for list command"""
    if string == 'all':
        return ('all', 0)
    elif string.isdigit():
        return ('num', int(string))
    elif string[-1] == '%':
        if string[:-1].isdigit():
            value = int(string[:-1])
            if value >= 0 and value <= 100:
                return ('%', value)
    msg = "parameter value must be n, n% or 'all' where n is a valid integer"
    raise argparse.ArgumentTypeError(msg)

def error_exit(msg):
    """print error message in same format as argparse then exit script"""
    print(sys.argv[0] + ': error: ', msg, sep='')
    exit()

if __name__ == '__main__':
    main()
