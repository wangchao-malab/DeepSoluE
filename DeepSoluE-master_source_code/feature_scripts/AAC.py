#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
from collections import Counter

def AAC(fastas, **kw):
    myOrder = 'ACDEFGHIKLMNPQRSTVWY'
    kw = {'order': myOrder}
    AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    #AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
                   count[key] = count[key]/len(sequence)
        code = [name]
        for aa in AA:
                  code.append(count[aa])
        encodings.append(code)
    return encodings
def AAC_feature(fastas):
    import sequence_read_save
    encodings = AAC(fastas)
    sequence_read_save.save_to_csv(encodings, "AAC.csv")