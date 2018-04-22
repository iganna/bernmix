import numpy as np
from urllib.request import urlretrieve
from urllib.request import urlopen



def kgml_downloader(path_to_download):

    for pathway_hsa_id in range(10000):

        url = 'http://rest.kegg.jp/get/hsa' +\
              '{0:05d}'.format(pathway_hsa_id) +'/kgml'
        try:
            urlopen(url)
        except:
            continue
        print(pathway_hsa_id)
        filename = path_to_download + 'hsa{0:05d}'.format(pathway_hsa_id) + '.xml'
        urlretrieve(url, filename)



if __name__ == "__main__":
    path_to_download = 'data/kgml/'
    kgml_downloader(path_to_download)

