#!/usr/bin/python

import sys
import os
import time
import re
import synapseclient
from synapseclient import File
import logging
log = logging.getLogger(__name__)
out_hdlr = logging.StreamHandler(sys.stdout)
out_hdlr.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
out_hdlr.setLevel(logging.INFO)
log.addHandler(out_hdlr)
log.setLevel(logging.INFO)

syn = synapseclient.login()

tbl = syn.tableQuery("SELECT specimenID,individualID,name,id,assay,dataType,sex,dataSubtype,consortium,study,diagnosis,tumorType,isMultiIndividual,isMultiSpecimen,isCellLine,species,fundingAgency,resourceType,nf1Genotype,nf2Genotype,studyName FROM syn11614202 WHERE ( ( \"assay\" = 'rnaSeq' ) AND ( \"fileFormat\" = 'fastq' ) AND ( \"parentId\" = 'syn7979306' ) ) order by specimenID")

df = tbl.asDataFrame()

#format = lambda x: x.replace('.fastq.gz','').split('_')[5]
#df['samp'] = df['name'].map(format)

#first run the index file
gencodeV29=syn.get('syn18134565').path
while not os.path.exists(gencodeV29):
    time.sleep(1)
ind_cmd='./salmon index --gencode -t '+gencodeV29+' -i gencode_v29_index'
print(ind_cmd)
if not os.path.exists('gencode_v29_index'):
    os.system(ind_cmd)


grouped = df.groupby('specimenID')
for name, group in grouped:
    file_list_f1 = []
    file_list_f2 = []
    for index,row in group.iterrows():
        file_path = syn.get(row['id']).path

        while not os.path.exists(file_path):
            time.sleep(1)
        if re.search('_1_GSL',file_path):##this needs to be customized, ideally
                                        ##add mate 1/2 to annotations
            file_list_f1.append(file_path)
        else:
            file_list_f2.append(file_path)

    # run salmon
    scmd='./salmon quant -i gencode_v29_index -l A -1 '+" ".join(file_list_f1)+' -2 '+" ".join(file_list_f2)+' -o '+os.path.join("quants",name)
    print(scmd)
    log.debug(scmd)
    if not os.path.exists(os.path.join('quants',name,'quant.sf')):
        os.system(scmd)


# upload to Synpase
df.drop(columns=['id', 'name'],inplace=True)
df.drop_duplicates(inplace=True)
df.reset_index(drop=True,inplace=True)
df['dataSubtype'] = 'processed'
df['fileFormat'] = 'sf'

format = lambda x: x.replace(' ','_')
df = df.fillna('')

for index,row in df.iterrows():
    folder_name = format(row['specimenID'])
    annotations = row.to_dict()
    temp = File(path=os.path.join("quants",folder_name,"quant.sf"), name="_".join([folder_name,"Salmon_gencodeV29","quant.sf"]), annotations=annotations,parent='syn18407530')
    temp = syn.store(temp)
    print(temp.id)
