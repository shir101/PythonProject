from DataCollector import Covid19DataPortalAccessionFetcher as cdc
import shutil
import os

#dc = cdc()  # instance of the data collector
print("dc cdc successful!!!!!")
#avail_lineages = dc.getLineagesList()  ####
avail_lineages = ['B.1.1', 'B.1.1.285', 'B.12', 'B.1.1.284', 'B.1.1.214']
# print(avail_lineages)
source = "./data/raw/"
dest = "./data/"

for lineage in avail_lineages[1:3]:
    try:
        os.mkdir("./data/"+lineage)
        dest2 = "./data/"+lineage
        # accessions = dc.getAccessionsByLineage(lineage)  #### take accessions of a certain lineage
        # dc.downloadAccessions(accessions)
        files = os.listdir(source)
        for file in files:
            shutil.move(os.path.join(source,file),dest2)
    except FileExistsError:
        pass