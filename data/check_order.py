#!/home/jrca253/anaconda2/bin/python

fTrain = open('tmp1.txt', 'r')
fValid = open('tmp2.txt', 'r')

train_cpgs = fTrain.readlines()
valid_cpgs = fValid.readlines()

fTrain.close()
fValid.close()

order_matches = True

for i in range(0, len(train_cpgs)):
    if i % 100000 == 0:
        print "Processing line " + str(i+1)
    cpg_t = train_cpgs[i].replace('\n', '')
    cpg_v = valid_cpgs[i].replace('\n', '')

    if cpg_t == cpg_v:
        continue
    else:
        print "MISMATCH at line " + str(i+1)
        order_matches = False
        break
    #ENDIF

#ENDFOR

if order_matches:
    print "CpG order matches in both files!"
