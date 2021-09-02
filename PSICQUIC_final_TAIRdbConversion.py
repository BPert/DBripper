from operator import index
import requests  # python -m pip install requests
import os
import subprocess
from subprocess import Popen,PIPE
import argparse
from collections import defaultdict
from more_itertools import unique_everseen
import time
from time import sleep
from progress.bar import Bar #pip install progress
#from progress.spinner import Spinner #pip install progress
#from progress.bar import ProgressBar
import shutil
#from bioservices import PSICQUIC
import pandas as pd
import openpyxl
from requests.models import HTTPError
import xlsxwriter
import re
# from tqdm import tqdm
from datetime import datetime
import sys
import urllib.parse
import urllib.request
import json


########################  Creating a parser

parser = argparse.ArgumentParser(prog='DBripper',
                                 usage='%(prog)s [-h] [-v] [-s] [-vurl] [-allmt] [-sp] [-dm] [-ch] [-it] [-pubselect] [-Prot] prot1 prot2... OR/AND @filename.txt ',
                                 description='Process some proteins. For -dm, -it and -ch use the MI ontology without the "MI:", just the four number.',
                                 fromfile_prefix_chars='@')

                                
parser.add_argument('-Prot',
                    metavar='prot',
                    help="enter all your proteins of interest, separated by a space OR/AND with a .txt file with @filename.txt" ,
                    nargs='*')
                    #required=True) #c'est pas valid pour pour les positionnels!
parser.add_argument('-dm',
                    '--detection_method',
                    type=str,
                    nargs='+',
                    dest='dm',
                    action='store',
                    default=None,
                    help= 'PSI-MI Interaction Detection Method, ie(two hybrid ): -dm 0018 ')
parser.add_argument('-it',
                    '--interaction_type',
                    type=str,
                    nargs='+',
                    dest='it',
                    action='store',
                    default=None,
                    help= 'PSI-MI Interaction Type, ie(physical association ): -it 0915 ')
parser.add_argument('-sp',
                    '--species',
                    type=str,
                    nargs='+',
                    dest='sp',
                    action='store',
                    default=str(3702),
                    help= 'The species you want to search in. Default : 3702(Arabidopsis thaliana)')
parser.add_argument('-db',
                    '--database',
                    type=str,
                    nargs='+',
                    dest='db',
                    action='store',
                    default=None,
                    help= 'Select the database to search in. First you can execute the script without -db and see all the correct name of the actives DBs in the psicquicActiveDB.txt in the SAVED_FILES directory!')
parser.add_argument('-ch',
                    '--child',
                    type=str,
                    nargs='+',
                    dest='ch',
                    action='store',
                    default=None,
                    #const="all",
                    help= 'Add to your detection methods selected all the descendants in the MI ontology, ie(two hybrid): -ch 0018 .\n If 0018 is selected in -dm it will also filter your results by all the descendants of 0018. If you want all descendants from all detection methods selected : -ch all ')
parser.add_argument('-pubselect',
                    '--publication_selection',
                    action='store_true',
                    help='save only one interaction for the same pubmed ID.')
parser.add_argument('-v',
                    '--verbose',
                    action='store_true',
                    help='an optional argument')
parser.add_argument('-vurl',
                    '--verbose_url',
                    action='store_true',
                    help='an optional argument')
parser.add_argument('-s',
                    '--save',
                    action='store_true',
                    help='Save all the created files in the "SAVE_FILES" directory')
parser.add_argument('-allmt',
                    '--allmethodesandtypes',
                    action='store_true',
                    help='Use this argument to know all the detection methods and all the interaction types encountered for your queried protein')

args = parser.parse_args()


print("args.Prot c'est ",args.Prot)



'''##########################  aaaaaaaaaaaaaaaaaaaaaaa    VVVVVVVIIIIRRREER
lgenes=[]
lgenes=args.Prot
print("la lgenes est",lgenes)

if args.verbose:
    print("\nYour list of query proteins:"," ".join(args.Prot),"\n")
sys.exit()
##########################  aaaaaaaaaaaaaaaaaaaaaaa    VVVVVVVIIIIRRREER'''










########################## Create a new args.dm if the argument -ch (--child) is selected. 
#                          It will contain all the descendants (in the MI ontology) of the detection methodes selected 

print("le args.dm est :",args.dm)
print("le args.ch est :",args.ch)

if args.ch != None:
    if args.ch[0] == "all":
        parentList=args.dm
    if args.ch and args.ch[0] != "all":
        parentList=args.ch

    print("\n\n LA parentList est:",parentList,"\nle type(parentList)",type(parentList))

dicoChild=defaultdict(str)

if not os.path.exists('JSON_OLS'):
    os.makedirs('JSON_OLS')

skip=False
if args.ch != None:
    try:
        for parent in parentList:
            url='http://www.ebi.ac.uk/ols/api/ontologies/mi/hierarchicalDescendants?id=MI_'+parent+'&size=500'

            print("\n\n l'url est :",url)
            req_response = requests.get(url)

            with open ("JSON_OLS/reponse_MI.json","wb") as f1:
                f1.write(req_response.content)
            
            with open("JSON_OLS/reponse_MI.json","r") as a:
                try:
                    dict1 = json.load(a)
                except json.decoder.JSONDecodeError:
                    print("It seems that your -ch MI number is incorrect or not found. Retry the command line.")
            try:
                for i in range(len(dict1["_embedded"]["terms"])):
                    dicoChild[dict1["_embedded"]["terms"][i]["annotation"]["id"][0].split("I:")[1]]=dict1["_embedded"]["terms"][i]["label"]
            except KeyError:
                print("\n The detection method ",parent," as no descendant.")
    except TypeError:
        print("The parameter -ch is misused.")
        skip=True
    try:
        args.dm = [*args.dm,*dicoChild.keys()]
    except TypeError:
        if skip==False:
            print("The parameter -ch is misused.")
        sys.exit()

######################## Class creat??????????

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '3[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

######################## STRING  searching string_ID

dicoProtID = defaultdict(str) # dictionnary keys = input id protein  /  value = the corresponding string ID

if args.verbose:
    print(f"{color.BOLD}{color.GREEN}{color.UNDERLINE}Searching the String Database.{color.END}\n")

# HTTPS request

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"
request_url = "/".join([string_api_url, output_format, method])

for prot in args.Prot:
    # Set parameters

    params = {

        # your protein list
        "identifiers": prot, # the input prot
        "species": args.sp,  #  organism identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
    }

    # Query STRING

    results = requests.post(request_url, data=params)

    # Record the results in a dictionary dicoProtID
    try:
        for line in results.text.strip().split("\n"):
            l = line.split("\t")
            dicoProtID[l[0]]=l[2]
            input_identifier, string_identifier = l[0], l[2]
        '''if args.verbose:
            print("Input:", input_identifier, "STRING_ID:", string_identifier, sep="\t")'''
    except IndexError:
        dicoProtID[prot]=prot
        if args.verbose:
            print(f"{color.BOLD}{color.RED}{color.UNDERLINE}PAS de correspondance string ID pour: {prot}\n{color.END}")

for prot in dicoProtID:
    if args.verbose:
        if dicoProtID[prot] == prot:
            print("Input:",prot,"value is also the input:",dicoProtID[prot],"string ID not found!")
        else:
            print("Input:", prot, "string ID:", dicoProtID[prot], sep="\t")
if args.verbose:
    print()



################################# STRING requests and retrieving DATA...

string_api_url = "https://string-db.org/api"
output_format = "psi-mi-tab"
method = "network"

# Construction de l' URL

request_url = "/".join([string_api_url, output_format, method])

# Pour chaque protéine une requête à STRING

for prot in dicoProtID:

    # Set parameters
    
    params = {

        "identifiers": dicoProtID[prot],  # String ID protein OR prot input if no string ID found.
        "species": args.sp,  # species  identifier
        #"caller_identity": "www.awesome_app.org"  # your app name
    }
    
    # Call STRING

    ###  requêtes via PSICQUIC  :  request_url = "http://string.uzh.ch/psicquic/webservices/current/search/interactor/"+dicoProtID[prot]+"?format=tab25"
    ### si utilisée >> enlever construction url et params!!!!!!!
    
    response = requests.post(request_url,data=params)

    # Save the network to file

    file_name = "%s_string_network_API.txt" % prot
    if args.verbose_url:
        print("Saving interaction network to %s" % file_name)

    if not os.path.exists('STRING_Data'):
        os.makedirs('STRING_Data')

    with open("./STRING_Data/"+file_name, 'wb') as fh:
        fh.write(response.content)

    file_name_out = file_name.split('.')[0]+'_sansDoubles.txt'
    if args.verbose_url:
        print("and saving to",file_name_out)
    with open("./STRING_Data/"+file_name, 'r') as f, open("./STRING_Data/"+file_name_out, 'w') as out_file:
        out_file.writelines(unique_everseen(f))
    
    if args.verbose_url: # processing bar  -> from progress.bar import Bar
        with Bar('Processing', max=4) as bar:
            for i in range(4):
                time.sleep(0.25)
                bar.next()
            print()

            
#######################

sp = subprocess.run(['Rscript','rassembleString_api.R'],stderr=subprocess.PIPE)

###########################

if args.verbose:
    print(f"{color.BOLD}{color.GREEN}{color.UNDERLINE}Searching with PSICQUIC.{color.END}\n")
    print(f"{color.UNDERLINE} Searching for UniProtKB ID with UniProt API{color.END}\n")
    
########################### PSICQUIC - searching for uniproKB ID using UniProt API with STRING_ID


dicoProtUniID=defaultdict(str) # dictionnary keys = input id protein  /  value = the corresponding uniprotKb ID
lnotFoundUniprot=[]
for prot in dicoProtID: # dictionnary keys = input id protein  /  value = the corresponding string ID
    
    url = 'https://www.uniprot.org/uploadlists/'
    
    paramsuniID = {
    'from': 'STRING_ID', # the input -> abbreviation of the ID format
    'to': 'ACC', # the output -> abbreviation of the ID format (uniprotKB)
    'format': 'tab', # for the output file
    'query': dicoProtID[prot]
    }

    data = urllib.parse.urlencode(paramsuniID)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    
    try:
        with urllib.request.urlopen(req) as f:
            responseuniID = f.read()
    #print('la prot en input',prot,'la query',dicoProtID[prot] ,' la réponse',responseuniID.decode('utf-8'))
    #print(" la reponse.split('tab')[2]",responseuniID.decode('utf-8').split('\t')[2].rstrip('\n'))
            try:
                
                dicoProtID[prot]=responseuniID.decode('utf-8').split('\t')[2].rstrip('\n')
                #print('coucoi le dicoProtUniID[prot] est ',dicoProtUniID[prot])
            except IndexError:
                dicoProtID[prot]=prot
                lnotFoundUniprot.append(prot)
                if args.verbose:
                    print(f"{color.BOLD}{color.RED}{color.UNDERLINE}No UniProtKb ID found for: {dicoProtID[prot]}\n{color.END}")

    except OSError as err:
        dicoProtID[prot]=prot
        if args.verbose:
            print(f"{color.BOLD}{color.RED}{color.UNDERLINE}> ERROR  while searching for uniprotkb ID with: {dicoProtID[prot]} / input protein:{prot}{color.END}")
            print("OS error: {0}".format(err),"\n")

for prot in dicoProtID:
    if args.verbose:
        if dicoProtID[prot] == prot:
            print("Input:",prot,"query request",dicoProtID[prot],"value is also the input:",dicoProtID[prot],"uniprotkb ID not found!")
        else:
            print("Input:", prot,"query request",dicoProtID[prot], "uniprotKB:", dicoProtID[prot], sep="\t")
if args.verbose:
    print()

########################### PSICQUIC - searching for uniproKB ID using TAIR conversion file

if args.verbose:
    print(f"{color.UNDERLINE} Searching for UniProtKB ID with TAIR file{color.END}\n")

lnotFoundTAIR=[]
lgenes=args.Prot
lgfound=[]
dicoProtIdTairF=defaultdict(str)

with open ('TAIR2UniprotMapping.txt','r') as f1:
	for l in f1:
		lp=l.rstrip('\n')
		lps=lp.split('\t')
		if lps[2] in lgenes:
			lgfound.append(lps[2])
			print (" Input TAIR - ",lps[2]," -> UniProtKB - ",lps[0])
			dicoProtIdTairF[lps[2]]=lps[0]

for gene in lgenes:
    if gene not in lgfound:
        if args.verbose:
            print(f"{color.BOLD}{color.RED}{color.UNDERLINE}No UniProtKb ID found for: {gene}\n{color.END}")
        lnotFoundTAIR.append(gene)

######################### Keeping the best way to convert TAIR ID into UniprotKb ID

if len(lnotFoundTAIR) < len(lnotFoundUniprot):
    if args.verbose:
        print("More UniProtKb found with the TAIR file")
    dicoProtID=dicoProtIdTairF

if len(lnotFoundTAIR) >= len(lnotFoundUniprot):
    if args.verbose:
        print("More UniProtKb found with the UniProt API using STTRING_ID")





'''    print("  333333333333333 changement de dico")

print("   33333333333  le dicoProtIdTAIRF est ",dicoProtIdTairF)
print(" 333333333333333  le dicoProt est ",dicoProtID)

sys.exit()'''


########################### PSICQUIC searching for ACTIVE DB...

dicoActiveDB=defaultdict(str) # dictionnary key : active Database name / value : url

psicquicActiveDB_response = requests.get(
    "http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=txt&protocol=rest"
)

now = datetime.now()
dt_string = now.strftime("%d.%m.%Y_%H:%M:%S")


with open ('psicquicActiveDB_'+dt_string+'.txt','wb') as f1:
    f1.write(psicquicActiveDB_response.content)

with open ('psicquicActiveDB_'+dt_string+'.txt','r')as f1:
    if args.verbose:
        print('All the active Databases:')
    for line in f1:
        lp=line.rstrip('\n').split('=')
        dicoActiveDB[lp[0]]=lp[1]
        if args.verbose:
            print(lp[0], end=" ")

removedb=[]

if args.db:
    for db in args.db:
        if db not in dicoActiveDB.keys():
            print(f'{color.RED}\n\n The database - {db} - selected is not active or simply misspelled (see above "All the active Databases") ')
            sys.exit()
    for db in dicoActiveDB:
        if db not in args.db:
            removedb.append(db)
    for db in removedb:
        del dicoActiveDB[db]
    print(f'\n\n The selected Database :{" ".join(dicoActiveDB.keys())}')



########################### PSICQUIC requests on all the actives DB found and retrieving DATA ...

if args.verbose :
    print("\n\nPSICQUIC requests...\n\n")
   
listeFile=[]

for activeDB in dicoActiveDB:
    for prot in  dicoProtID:

        try:
            reqACDATA_response = requests.get(
                dicoActiveDB[activeDB]+'query/id:'+dicoProtID[prot]+' AND species:'+args.sp+'?format=tab25')
            if args.verbose_url:
                print(f'DB :{color.BLUE}{activeDB}{color.END}\nREQUEST URL={dicoActiveDB[activeDB]}query/id:{dicoProtID[prot]} AND species:{args.sp}?format=tab25')
                
        except IOError:
            if args.verbose:
                print(f'{color.RED}DB :{activeDB}\nCan not open url:{dicoActiveDB[activeDB]}query/id:{dicoProtID[prot]} AND species:{args.sp}?format=tab25{color.END}')

        
        file_name = activeDB+"__"+prot+".txt"

        listeFile.append(file_name)

        if not os.path.exists('PSICQUIC_Data'):
            os.makedirs('PSICQUIC_Data')
        with open ("./PSICQUIC_Data/"+file_name,'wb') as f1:
            f1.write(reqACDATA_response.content)
        with open ("./PSICQUIC_Data/"+file_name,'r') as f1:  # to add the column with the TAIR ID.
            lines=f1.read()
            lines=lines.replace('\n','\t'+prot+'\n')
        with open ("./PSICQUIC_Data/"+file_name,'w') as f1:
            f1.write(lines)


################################# Trimming wrong files (empty ones, error messages from databases)

# Trim all the empty files ; no data given by the request.

dicoNbFileProt=defaultdict(list) # dictionnary keys = protein / value : the liste of files used (with interaction found) for this protein
                                #   if no interaction found the list is empty
listeNotEmpty=[]
for fichier in listeFile:
    protid=fichier.split('__')[1].split('.')[0]
    taille=os.path.getsize("./PSICQUIC_Data/"+fichier)
    if taille != 0:
        listeNotEmpty.append(fichier)
    else:
        dicoNbFileProt[protid]

# Filter all the files that are not in MITAB 2.5 format ; those corresponding to error message from data base.
newListFile=[]

if args.verbose:
    print(f'{color.BOLD}{color.RED}\n UNUSED FILES (because of error message from the data base or unknown format):\n{color.END}')
for fichier in listeNotEmpty:
    protid=fichier.split('__')[1].split('.')[0]
    with open ("./PSICQUIC_Data/"+fichier,'r') as f1:
        ctlok=0
        ctl=0
        for l in f1:  # to check if the file is MITAB 2.5 format
            ctl+=1
            lp=l.rstrip('\n').split('\t')
            if len(lp) == 16:
                ctlok+=1
        if ctlok == ctl:
            newListFile.append("./PSICQUIC_Data/"+fichier)
            dicoNbFileProt[protid].append(fichier)
        else:
            dicoNbFileProt[protid]
            if args.verbose:
                print(f"{color.RED}{fichier}{color.END}")
if args.verbose:
    if len(newListFile) == len(listeNotEmpty):
        print(f"{color.RED} -> None. All files are used.{color.END}")

####################### Write the file TOTAL_Psicquic_nonFILTRE


with open('./PSICQUIC_Data/TOTAL_Psicquic_nonFILTRE.txt','wb') as wfd:
    for f in newListFile:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd, 1024*1024*10)

tailleTXT=os.path.getsize('./PSICQUIC_Data/TOTAL_Psicquic_nonFILTRE.txt')

if tailleTXT == 0:
    print(f'{color.RED}{color.BOLD}\n\nNO INTERACTION FOUND WITH ANY OF YOUR QUERY PROTEIN:{" ".join(args.Prot)}\n{color.END}')
    sys.exit()


df = pd.read_table('./PSICQUIC_Data/TOTAL_Psicquic_nonFILTRE.txt', header=None)

df.to_excel('./PSICQUIC_Data/TOTAL_Psicquic_nonFILTRE.xlsx', 'Sheet1', index=False, header=False)

######################## FILTRER les doublons!

dicoTOT=defaultdict(lambda:defaultdict(str))
newMeth=[] #a list that will contain all the different detection methods encountered in all the dataset IF NO ARGS dm/it
newType=[] #a list that will contain all the different interaction types encountered in all the dataset IF NO ARGS dm/it

def New(numlp,param):  # a fonction to create a list of detection methods OR a list of interaction types encountered IF NO ARGS dm/it
    numIM=re.search('MI:[0-9]{4}',lp[numlp]).group().split(':')[1]+":"+re.search('\(.*\)',lp[numlp]).group()
    param.append(numIM) if numIM not in param else param

# The goal of  following 'with open' is to create a dictionnary of dictionnary (dicoTOT) 
#   key : protein ID 1 +"\t"+ prot ID2 +"_/_n"+the number of time that the protein couple has been encountered
#   value : three new key "meth" "typ" "restcol":
#                -key "meth": the detection method for this couple only
#                -key "typ": the interaction type for this couple only
#                -key "restcol": a list containing all the columns except col A and col B 
ctDoublons=0
dicoDJV=defaultdict(int)
with open("./PSICQUIC_Data/TOTAL_Psicquic_nonFILTRE.txt",'r') as f1txt:
    for l in f1txt:
        lp=l.rstrip("\n").split("\t")
        if lp[0]+"\t"+lp[1]+"_/_n1" in dicoTOT.keys():
            ctDoublons+=1
            if lp[0]+"\t"+lp[1] not in dicoDJV.keys():
                '''#print("\n\nON PASSE  PAR LE NONE§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\n\n")'''
                dicoDJV[lp[0]+"\t"+lp[1]]=1
            dicoDJV[lp[0]+"\t"+lp[1]]+=1
            dicoTOT[lp[0]+"\t"+lp[1]+"_/_n"+str(dicoDJV[lp[0]+"\t"+lp[1]])]["restcol"]=lp[2:]

            if lp[6] != "-":
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n"+str(dicoDJV[lp[0]+"\t"+lp[1]])]["meth"]=re.search('MI:[0-9]{4}',lp[6]).group().split(':')[1]
                New(6,newMeth)
            else:
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n"+str(dicoDJV[lp[0]+"\t"+lp[1]])]["meth"]="-"
            if lp[11] != "-":
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n"+str(dicoDJV[lp[0]+"\t"+lp[1]])]["typ"]=re.search('MI:[0-9]{4}',lp[11]).group().split(':')[1]
                New(11,newType)
            else:
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n"+str(dicoDJV[lp[0]+"\t"+lp[1]])]["typ"]="-"
        else:
            dicoTOT[lp[0]+"\t"+lp[1]+"_/_n1"]["restcol"]=lp[2:]

            if lp[6] != "-":
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n1"]["meth"]=re.search('MI:[0-9]{4}',lp[6]).group().split(':')[1]
                New(6,newMeth)
            else:
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n1"]["meth"]="-"
            if lp[11] != "-":
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n1"]["typ"]=re.search('MI:[0-9]{4}',lp[11]).group().split(':')[1]
                New(11,newType)
            else:
                dicoTOT[lp[0]+"\t"+lp[1]+"_/_n1"]["typ"]="-"

listMethArgs=[]
listTypeArgs=[]

def MethTypeArgs(mORt,listMorIT): # IF ARGS: create list of DM or IT encountered in all the dataset IF ARGS!
    for ids in dicoTOTmethORit:
        listMorIT.append(dicoTOTmethORit[ids][mORt]) if dicoTOTmethORit[ids][mORt] not in listMorIT else listMorIT



dicoTOTmethORit=defaultdict(lambda:defaultdict(str)) # to create a dictionnary regarding the args selected:
#                                                        Detection Methods OR/AND Interaction Type 

if args.dm and not args.it:
    for ids in dicoTOT:
        if dicoTOT[ids]["meth"] in args.dm:
            dicoTOTmethORit[ids]["restcol"]=dicoTOT[ids]["restcol"]
            dicoTOTmethORit[ids]["meth"]=dicoTOT[ids]["meth"]
            dicoTOTmethORit[ids]["typ"]=dicoTOT[ids]["typ"]
            

if args.it and not args.dm:
    for ids in dicoTOT:
        if dicoTOT[ids]["typ"] in args.it:
            dicoTOTmethORit[ids]["restcol"]=dicoTOT[ids]["restcol"]
            dicoTOTmethORit[ids]["meth"]=dicoTOT[ids]["meth"]
            dicoTOTmethORit[ids]["typ"]=dicoTOT[ids]["typ"]



if args.dm and args.it:
    for ids in dicoTOT:
        if  dicoTOT[ids]["meth"] in args.dm and dicoTOT[ids]["typ"] in args.it:
            dicoTOTmethORit[ids]["restcol"]=dicoTOT[ids]["restcol"]
            dicoTOTmethORit[ids]["meth"]=dicoTOT[ids]["meth"]
            dicoTOTmethORit[ids]["typ"]=dicoTOT[ids]["typ"]
          

MethTypeArgs("meth",listMethArgs)
MethTypeArgs("typ",listTypeArgs)

if args.allmethodesandtypes :
    if args.dm == None and args.it == None:
        print(f'\n{color.BOLD}{color.UNDERLINE}With no args selected:{color.END}\n -All detection methods encountered are: \n{newMeth}\n -All interaction types encountered are:\n{newType}')
    if args.dm or args.it:
        print(f'\n{color.BOLD}{color.UNDERLINE}With this list of agrs selected:{color.END}\n- Detection methods:{args.dm}\n- Interaction types:{args.it}\n -All detection methods encountered are: \n{listMethArgs}\n -All interaction types encountered are:\n{listTypeArgs}')


dicoMethCounter=defaultdict(list)
dicoTypCounter=defaultdict(list)


def methCounter(**param): #from the selected dictionnary to create a dictionnary : keys = couple protein  
                            #/  values = liste of DM encountered for this couple protein
    for ids in param:
        idsp=ids.split("_/_n")[0]
        if param[ids]["meth"] not in dicoMethCounter[idsp] and param[ids]["meth"] != '-':
            dicoMethCounter[idsp].append(param[ids]["meth"])
        if param[ids]["typ"] not in dicoTypCounter[idsp] and param[ids]["typ"] != '-':
            dicoTypCounter[idsp].append(param[ids]["typ"])       

def dicoTOTcountPlus(**param): # to add on the selected dictionnary (in restcol) 
                                #the number of differents DM encountered for one protein and then the list of those different meth. encountered
    for ids in param:
        idsp=ids.split('_/_n')[0]
        if not dicoMethCounter[idsp]:
            param[ids]["restcol"].append('-')
            param[ids]["restcol"].append('-')
        else:
            param[ids]["restcol"].append(len(dicoMethCounter[idsp]))
            param[ids]["restcol"].append(dicoMethCounter[idsp])
        if  not dicoTypCounter[idsp]:
            param[ids]["restcol"].append('-')
            param[ids]["restcol"].append('-')
        else:
            param[ids]["restcol"].append(len(dicoTypCounter[idsp]))
            param[ids]["restcol"].append(dicoTypCounter[idsp])        
            

dicoTOTsansInv=defaultdict(lambda:defaultdict(str))
listKeysSansN=[i.split("_/_n")[0] for i in dicoTOT.keys()]  
dicoCountSameId=defaultdict(int)
for c in listKeysSansN:
    dicoCountSameId[c]=dicoCountSameId.get(c,0)+1

 
def sansInversion(**param):
    listIdDJV=[] 

    for ids in param:
        idspropres=ids.split('_/_')[0].split('\t')
        if idspropres[0] == idspropres[1]:   # to not reverse if prot A = prot B!!
            dicoTOTsansInv[ids]=param[ids]
            dicoTOTsansInv[ids]['restcol'].append('-')
        else:
            idSN=ids.split('_/_n')[0]
            idsinvSN=idSN.split("\t")[1]+"\t"+idSN.split("\t")[0]
            if ids not in listIdDJV :
                dicoTOTsansInv[ids]=param[ids]
                dicoTOTsansInv[ids]['restcol'].append('-')            
                listIdDJV.append(ids)
                for key in param:
                    if key.split('_/_n')[0] == idsinvSN : 
                        if key not in listIdDJV:
                            dicoCountSameId[idSN]+=1
                            idsinvrev=idSN+"_/_n"+str(dicoCountSameId[idSN])
                            dicoTOTsansInv[idsinvrev]=param[key]      
                            dicoTOTsansInv[idsinvrev]['restcol'].append('rev')
                            listIdDJV.append(key)



if args.dm or args.it:
    methCounter(**dicoTOTmethORit)
    dicoTOTcountPlus(**dicoTOTmethORit) 

if args.dm == None and args.it == None:
    methCounter(**dicoTOT)
    dicoTOTcountPlus(**dicoTOT)


if args.dm or args.it:
    sansInversion(**dicoTOTmethORit)
    
if args.dm == None and args.it == None:
    sansInversion(**dicoTOT)


############################ Convertion of DOI into PMID

if not os.path.exists('JSON_FILES'):
    os.makedirs('JSON_FILES')

for ids in dicoTOTsansInv:
    pubids=dicoTOTsansInv[ids]["restcol"][6].split("|")
    for pubid in pubids:
        if pubid.startswith('pubmed:'):
            dicoTOTsansInv[ids]["restcol"][6]=pubid
        if pubid.startswith('DOI:'):
            url="https://api.fatcat.wiki/v0/release/lookup?doi="+pubid.split(":")[1]
            req_response = requests.get(url)

            with open ("JSON_FILES/"+pubid.split(":")[1].split("/")[0]+".json","wb") as f1:
                f1.write(req_response.content)
            
            with open("JSON_FILES/"+pubid.split(":")[1].split("/")[0]+".json","r") as a:
                dict1 = json.load(a)
            dicoTOTsansInv[ids]["restcol"][6]="pubmed:"+dict1["ext_ids"]["pmid"]


##################### to keep just one interaction of all interaction from the same pubmed ID


dicoDJV2=defaultdict(lambda:defaultdict(list))  # a dictionnary keys : ids of the two proteins interacting  
                                                # / value : two new keys:- '-'(not reverse) value: the list of pubmed id encountered
                                                #                        - 'rev' (reverse)  value: the list of pubmed id encountered

dicoTOTsN=defaultdict(lambda:defaultdict(str))




if args.publication_selection:

    for ids in dicoTOTsansInv:
        idsSN=ids.split("_/_n")[0]
        
        if dicoTOTsansInv[ids]["restcol"][18] == '-':
            if dicoTOTsansInv[ids]["restcol"][6] not in dicoDJV2[idsSN]['-']: 
                dicoDJV2[idsSN]['-'].append(dicoTOTsansInv[ids]["restcol"][6])
        
        if dicoTOTsansInv[ids]["restcol"][18] == 'rev':
            if dicoTOTsansInv[ids]["restcol"][6] not in dicoDJV2[idsSN]['rev']:
                dicoDJV2[idsSN]['rev'].append(dicoTOTsansInv[ids]["restcol"][6])

    

    for ids in dicoTOTsansInv:
        idsSN = ids.split("_/")[0]
        if dicoTOTsansInv[ids]["restcol"][18] == '-':
            dicoTOTsansInv[ids]["restcol"].append(dicoDJV2[idsSN]['-'])
            dicoTOTsansInv[ids]["restcol"].append(len(dicoDJV2[idsSN]['-']))
            dicoTOTsN[idsSN]["restcol"]=dicoTOTsansInv[ids]["restcol"]
        if dicoTOTsansInv[ids]["restcol"][18] == 'rev':
            dicoTOTsansInv[ids]["restcol"].append(dicoDJV2[idsSN]['rev'])
            dicoTOTsansInv[ids]["restcol"].append(len(dicoDJV2[idsSN]['rev']))
            oldID=idsSN.split("\t")[1]+"\t"+idsSN.split("\t")[0]
            dicoTOTsN[oldID]["restcol"]=dicoTOTsansInv[ids]["restcol"]


################### write the TOTAL_PSICQUIC.xlsx


if not args.publication_selection:
    fout = "./PSICQUIC_Data/TOTAL_PSICQUIC_FR.txt"
    with open(fout, "w") as dataPF:
        for ids in dicoTOTsansInv:
            rest="\t".join(str(case) for case in dicoTOTsansInv[ids]["restcol"])
            #print("rest :",rest)
            dataPF.write(ids+"\t"+rest)
            dataPF.write("\n")

if args.publication_selection:
    fout = "./PSICQUIC_Data/TOTAL_PSICQUIC_FR.txt"
    with open(fout, "w") as dataPF:
        for ids in dicoTOTsN:
            rest="\t".join(str(case) for case in dicoTOTsN[ids]["restcol"])
            #print("rest :",rest)
            dataPF.write(ids+"\t"+rest)
            dataPF.write("\n")


df = pd.read_table('./PSICQUIC_Data/TOTAL_PSICQUIC_FR.txt', header=None)

if not os.path.exists('PSICQUIC_TOTAL_Data'):
    os.makedirs('PSICQUIC_TOTAL_Data')

# The headers of the TOTAL_PSICQUIC.xlsx:
if not args.publication_selection:
    headerList=['Unique identifier for interactor A',
                'Unique identifier for interactor B',
                'Alternative identifier for interactor A',
                'Alternative identifier for interactor B',
                'Aliases for A',
                'Aliases for B',
                'Detection methods',
                'First author',
                'Identifier of the publication',
                'NCBI Taxonomy identifier for interactor A',
                'NCBI Taxonomy identifier for interactor B',
                'Interaction types',
                'Source databases',
                'Interaction identifier(s)',
                'Confidence score',
                'TAIR',
                'nb DM',
                'list of DM',
                'nb IT',
                'list of IT',
                'reverse(rev) or not(-)']

if args.publication_selection:
    headerList=['Unique identifier for interactor A',
                'Unique identifier for interactor B',
                'Alternative identifier for interactor A',
                'Alternative identifier for interactor B',
                'Aliases for A',
                'Aliases for B',
                'Detection methods',
                'First author',
                'Identifier of the publication',
                'NCBI Taxonomy identifier for interactor A',
                'NCBI Taxonomy identifier for interactor B',
                'Interaction types',
                'Source databases',
                'Interaction identifier(s)',
                'Confidence score',
                'TAIR',
                'nb DM',
                'list of DM',
                'nb IT',
                'list of IT',
                'reverse(rev) or not(-)',
                'pub list',
                'nb of diff pub']

listColToHide=[2,3,4,5,7,8,9,10,13,14]


writer = pd.ExcelWriter('./PSICQUIC_TOTAL_Data/TOTAL_PSICQUIC_'+dt_string+'.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1', index=False , header=headerList)

workbook  = writer.book
worksheet = writer.sheets['Sheet1']

# Add a header format.
header_format = workbook.add_format({
    'bold': True,
    'text_wrap': True,
    'valign': 'top',
    'fg_color': '#D7E4BC',
    'border': 1})

for col_num, value in enumerate(headerList):
    worksheet.write(0, col_num , value, header_format)
for col in listColToHide:   
    workbook.get_worksheet_by_name('Sheet1').set_column(col,col, None, None, {'hidden': 1})

writer.save()

##########################
ListeProtNoInter=[]
if args.verbose:
    if not all(len(dicoNbFileProt[prot]) != 0 for prot in dicoNbFileProt):
        print(f'{color.RED}{color.BOLD}\nProteins with no interactions found! :{color.END}')  

    for prot in dicoNbFileProt:
        if len(dicoNbFileProt[prot]) == 0:
            ListeProtNoInter.append(prot)
            print(f'{color.RED} -> {prot}{color.END}')


################## Write a .txt file with the command line and all the detection methods and interaction types encountered in all the total dataset

with open("Resume_"+dt_string+".txt", "w")as filout:

    filout.write("Your list of query proteins:\n")
    filout.write(" ".join(args.Prot))
    filout.write("\n")
    if args.dm:
        filout.write("\nDetection methods selected.\n")
        for dm in args.dm:
            filout.write(dm)
            filout.write("\n")
    else:
        filout.write("\nNo Detection methods selected.\n")
    if args.it:
        filout.write("\nInteraction types selected:\n")
        for it in args.it:
            filout.write(it)
            filout.write("\n")
    else:
        filout.write("\nNO Interaction types selected.\n") 

    filout.write("\n\n The command line executed was: python3 ")
    filout.write(" ".join(sys.argv[:]))


    if args.dm == None and args.it == None:
        filout.write("\n\nWith no args selected:\n  -All detection methods encountered are:\n")
        for meth in newMeth:
            filout.write(meth)
            filout.write("\n")
        filout.write("\n  -All interaction types encountered are:\n")
        for typ in newType:
            filout.write(typ)
            filout.write("\n")
    if args.dm or args.it:
        filout.write("\n\nWith those agrs selected:\n")
        filout.write("\n\n  -All detection methods encountered are:\n")
        for meth in listMethArgs:
            filout.write(meth)
            filout.write("\n")
        filout.write("\n  -All interaction types encountered are:\n")
        for typ in listTypeArgs:
            filout.write(typ)
            filout.write("\n")
    
    filout.write('\n\nProtein(s) with no interaction found:\n')
    for prot in ListeProtNoInter:
        filout.write(prot)
        filout.write('\n')

#############################

workbook = xlsxwriter.Workbook("whichDataBase_"+dt_string+".xlsx")
worksheet = workbook.add_worksheet()
merge_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 2}) 

row = 0
col = 0
ctdeb=0 
ctfin=0 

for prot in dicoNbFileProt:
    if len(dicoNbFileProt[prot]) != 0:
        ctdeb=row
        for file in dicoNbFileProt[prot]:
            with open("./PSICQUIC_Data/"+file, 'r') as fd:
                nbofline = 0
                while fd.readline():
                    nbofline += 1
            worksheet.write(row, col,     prot, merge_format)
            worksheet.write(row, col + 1, file.split('__')[0])
            worksheet.write(row, col + 2, nbofline)
            ctfin=row
            row += 1

        if ctdeb != ctfin:
            worksheet.merge_range(ctdeb, 0, ctfin, 0, prot ,merge_format)

workbook.close()


###################


if args.save:

    #PSICQUIC - save all the files .txt (the response for one protein and for one database) in the SAVE_FILES directory and erase them from ./PSICQUIC_Data
    if not os.path.exists('SAVED_FILES/'+dt_string+'/PSICQUIC_Data'):
        os.makedirs('SAVED_FILES/'+dt_string+'/PSICQUIC_Data')
    source_dir = './PSICQUIC_Data' 
    target_dir = 'SAVED_FILES/'+dt_string+'/PSICQUIC_Data'
    file_names = os.listdir(source_dir)
    for file_name in file_names:
        shutil.move(os.path.join(source_dir, file_name), target_dir)

    # save all the files - TOTAL_PSICQUIC_date.xlsx - the liste of active databases - the resume
    if not os.path.exists('SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data'):
        os.makedirs('SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data')
    shutil.move('./PSICQUIC_TOTAL_Data/TOTAL_PSICQUIC_'+dt_string+'.xlsx', 'SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data')
    shutil.move('psicquicActiveDB_'+dt_string+'.txt', 'SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data')
    shutil.move("Resume_"+dt_string+".txt", 'SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data')
    shutil.move("whichDataBase_"+dt_string+".xlsx", 'SAVED_FILES/'+dt_string+'/PSICQUIC_TOTAL_Data')

    #STRING - save all the files .txt (the response for one protein and for one database) in the SAVE_FILES directory and erase them from ./STRING_Data
    if not os.path.exists('SAVED_FILES/'+dt_string+'/STRING_Data'):
        os.makedirs('SAVED_FILES/'+dt_string+'/STRING_Data')
    source_dir = './STRING_Data' 
    target_dir = 'SAVED_FILES/'+dt_string+'/STRING_Data'
    file_names = os.listdir(source_dir)
    for file_name in file_names:
        shutil.move(os.path.join(source_dir, file_name), target_dir)




print('\n\n')