# -*- coding: utf-8 -*-

import context

import seqanalyzer.sequence_lib as seqlib


def test_seqlib():
    print(seqlib.rc_seq('actg'))


def main():
    test_seqlib()


if __name__ == '__main__':
    main()