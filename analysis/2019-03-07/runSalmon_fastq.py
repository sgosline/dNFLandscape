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

tbl = syn.tableQuery("SELECT specimenID,individualID,name,id FROM syn11614202 WHERE ( ( \"assay\" = 'rnaSeq' ) AND ( \"fileFormat\" = 'fastq' ) AND ( \"parentId\" = 'syn7979306' ) ) order by specimenID")

df = tbl.asDataFrame()

format = lambda x: x.replace('.fastq.gz','').split('_')[5]
df['id'] = df['name'].map(format)

#first run the index file
gencodeV29=syn.get('syn18134565').path
while not os.path.exists(gencodeV29)
    time.sleep(1)
ind_cmd='salmon index --gencode -t '+gencodeV29+'--i gencode_v29_index'
os.system(ind_cmd)

grouped = df.groupby('id')
for name, group in grouped:
    file_list_f1 = []
    file_list_f2 = []
    for index,row in group.iterrows():
        file_path = syn.get(row['id']).path
#        file_path = row['name']
#        os.rename(temp.path, file_path)
        while not os.path.exists(file_path):
            time.sleep(1)
        if re.search('_1_GSL',file_path):
            file_list_f1.append(file_path)
        else:
            file_list_f2.append(file_path)
    # combine fastq files
#    f1=os.path.join(name+'_1.fastq.gz')
#    f2=os.path.join(name+'_2.fastq.gz')
#    os.system("cat "+" ".join(file_list_f1)+" > "+f1)
#    os.system("cat "+" ".join(file_list_f2)+" > "+f2)
#    while not os.path.exists(f1) or not os.path.exists(f2):
#        time.sleep(1)
    # remove fastq files

    # run salmon
    scmd='salmon quant -i gencode_v29_index -l A -1 '+" ".join(file_list_f1)+' -2 '+" ".join(file_list_f2)+' -o '+os.path.join("quants",name)
    log.debug(scmd)
    os.system(scmd)
    for f in file_list_f1:
        os.remove(f)
    for f in file_list_f2:
        os.remove(f)
        #   os.remove(f1)
 #   os.remove(f2)

# upload to Synpase
df.drop(columns=['id', 'name', 'dataFileHandleId'],inplace=True)
df.drop_duplicates(inplace=True)
df.reset_index(drop=True,inplace=True)
df['dataSubtype'] = 'processed'
df['fileFormat'] = 'sf'

format = lambda x: x.replace(' ','_')
df = df.fillna('')

for index,row in df.iterrows():
    folder_name = format(row['specimenID'])
    annotations = row.to_dict()
    temp = File(path=os.path.join("quants",folder_name,"quant.sf"), name="_".join([folder_name,"Salmon_gencodeV29","quant.sf"]), annotations=annotations,parent='syn17077846')
    temp = syn.store(temp)
    print(temp.id)
