import re
import numpy as np

with open('nuc_1-11_aligned.dbn', 'r') as f:
    for i, line in enumerate(f):
        if i % 3 == 0:
            tmp = []
        tmp.append(line)
        if i % 3 == 2:
            pname = re.sub(' ', '', tmp[0][2:5])

#            # Try to feed the gaps into jViz
#            tmp[1] = re.sub('-', 'A', tmp[1])
#            tmp[2] = ''.join([tmp[2][i] for i, p in enumerate(tmp[1]) if p != '-'])
#            tmp[1] = ''.join([p for i, p in enumerate(tmp[1]) if p != '-'])
            with open('nuc_'+pname+'.dbn', 'w') as ff:
                for lline in tmp:
                    ff.write(lline)
