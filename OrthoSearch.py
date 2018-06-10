# -*- coding: utf-8 -*-
"""
Created on Thu May 17 11:34:07 2018

@author: Ines
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Alphabet
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import Entrez
Entrez.email='ines.cf96@gmail.com'
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
import csv

#

#
def seq_NCBI(acNumber): #este coloca o genoma já num ficheiro fasta ATENÇAO: Tem de ser escrito sem espaços!!!
    acNumber= str(acNumber)
    acNumber=acNumber.split(',') #caso receba mais que um acNumber!
    
#    file=open('genomes.txt', 'w') #abrir um ficheiro para ser a database
#    file.close()
    
    dic={}
    for x in acNumber: #ir a cada acNumber
        print('Fetching data for', x, 'phage genome...', end='') # para o Done aparecer à frente dos 3 pontinhos 
        
        record = Entrez.read(Entrez.elink(dbfrom="nucleotide", id=x,db="protein"))
        # print(record)
    #    for linksetdb in record[0]["LinkSetDb"]:
    #      print(linksetdb["DbTo"], linksetdb["LinkName"], (linksetdb["Link"])) 

        seq=[]
        name=[]
        for link in record[0]["LinkSetDb"][0]["Link"]:
    #         print(link["Id"])
             handle=Entrez.efetch(db='protein', rettype='fasta', retmode='text', id=link["Id"])
             seq_record=SeqIO.read(handle, 'fasta')
             handle.close()
    #         print(seq_record)
             seq.append(str(seq_record.seq)) #colocar sequencia numa lista
             name.append(seq_record.id) #colocar id da proteina numa lista
        dic[x] = name
        print('Done')
      
###- Isto seria para sacar a sequencia no formato de nucleotido( forma mais simples)         
###    handle=Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=acNumber)
###    seq_record=SeqIO.read(handle, 'fasta')
###    handle.close()
###    seq=[]
###    name=[]
##
###- Isto seria para obter a sequencia nucleotidica e traduzi-la   
###print(seq_record.id, 'sequencia', seq_record.seq)
### translated = seq_record.seq.translate() 
###   seq.append(str(translated)) #colocar sequencia numa lista
###   name.append(seq_record.id) #colocar id do genoma numa lista
##    
#    
        f = open("genomes.txt", "a") #adicionar a sequencia e o nome ao ficheiro já existente
        for i in range(len(seq)):
            f.write(">|" + x + "|" + name[i] + "\n" + seq[i] + "\n")
        f.close()
        
#    
#    
    return dic
##
##------------------------------------------------------------------------------------
#
##Importante!!
#
#
#
def blast():
###MAKEBLASTDB --> CREATE DATABASE   
##    
    cmd= "makeblastdb -in genomes.txt -dbtype nucl"
    
    try:
        subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
        output = subprocess.check_output(cmd, shell=True)
        print(output.decode("iso-8859-1")) #para conseguir ver que a database foi criada
    
    except subprocess.CalledProcessError as e:
        raise RuntimeError("Comando '{}' \nErro (código {}): {}".format(e.cmd, e.returncode, e.output))

#    
#    
###RUN BLAST
##    
    blastp_cline = NcbiblastpCommandline(query='genomes.txt', subject='genomes.txt', evalue=0.001, outfmt= 5, out='genomes.xml')
    blastp_cline
#
##---------------nao interessa---------!   
##    NcbiblastxCommandline(cmd='blastx', out='testing_sequences.xml', outfmt=5, query='sequences.txt', db='nr', evalue=0.001)
##---------------nao interessa---------!  
#
##-------Isto serve para ver se o meu BLAST foi feito----:     
#    #print(blastp_cline)
##-------Isto serve para ver se o meu BLAST foi feito----:
    stdout, stderr = blastp_cline()


#PARSING BLAST OUTPUT

    result_handle = open("genomes.xml")
    blast_records = NCBIXML.parse(result_handle) #multiplas queries
    E_VALUE_THRESH = 0.004
    
    
    align_res = {}
    for blast_record in blast_records:
        
        blast_record_query = blast_record.query.split('|')
#        print('\nQuery Sequence', blast_record_query) 
        
        results=[]
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hsp_value= hsp.identities/hsp.align_length
#                print(hsp.identities/hsp.align_length ) #>95% daquelas que sao verdadeiramente iguais, só aceito as que tiverem isto assim
                
                if hsp.expect < E_VALUE_THRESH:
                    #title_alignment= str(alignment.title) #originalmente
                    title_alignment = alignment.title.split('|')
#                    print('****Alignment****')
#                    print('sequence:', title_alignment)
#                    print('length:', alignment.length)
#                    print('e value:', hsp.expect)
#                    print('score:', hsp.score)
#                    
                    #results.append(alignment.title[2].rstrip()) #era o original!!
                    results.append((title_alignment[2].rstrip(), hsp_value, None))
        
            
        if blast_record_query[2] not in align_res:
            align_res[blast_record_query[2]]= results
#    
#    print (align_res)
#    
#   
#    #Conferir se alinham realmente (se são mesmo iguais):
    clusters={} #chave id do cluster..valores-->proteinas desse cluster
    iD=0
    for key in align_res:
        if len(align_res[key]) > 1: #se tiver mais do que uma proteina como valor (ou seja, se tiver alinhamentos para além de ela propria)
            a= list(align_res[key]) #lista com valores da key em especifico
            for value in range(len(align_res[key])):
                if align_res[key][value][0] != key :#para descartar o alinhamento com ela propria
                    if align_res[key][value][1] < 0.95: #descartar tudo o que tem menos de 0.95
                        a.remove(align_res[key][value])#vou aos valores da key e removo essa proteina com a qual alinhou..visto que nao sao a mesma
                        continue
                    
                    values = [x[0] for x in align_res[align_res[key][value][0]]] #cria uma lista de valores de cada proteina que está em value (so o nome das proteinas)
                    if key not in values: #este ultimo ponto dá os valores do valor quando ele é chave! Entao, se a chave que lhe deu origem nao está nos seus valores, nao sao a mesma proteina
                        a.remove(align_res[key][value])
                        
            align_res[key]= a #atualizar dicionario só com a proteina e as que sao iguais a ela
            
            var= False #INICIALMENTE NÃO ESTÁ EM NENHUM CLUSTER 
            for cluster in clusters:
              if key in clusters[cluster]:
                  var=True #SE ESTIVER PASSA A TRUE
                  break
            if not var: #a var é falso --> é porque não está em nenhum cluster
                
                if len(align_res[key]) == 1:
                    clusters[key]=[key]
                    align_res[key][0] = (align_res[key][0][0], align_res[key][0][1], key)
                
                else:
                    ident= 'Cluster' + str(iD)
                    clusters[ident] = [x[0] for x in align_res[key]]
                    iD += 1
                    for value in range(len(align_res[key])):
                        align_res[key][value] = (align_res[key][value][0],align_res[key][value][1], ident)
                        
                        
                    
                
            
    #print (align_res) #este dicionario já é o que me interessa...aquele que possui a proteina e as iguais a ela!
    #print(clusters)
    return clusters


##CREATE MATRIX
def matrix(dic_trackback, dic_clusters):
    
    #Passo1: criar as listas
    genomes=[] #lista com todos os genomas
    proteins= list(dic_clusters.keys()) #lista com todas as proteinas
    
    for key in dic_trackback: 
        genomes.append(key) #colocar na lista de genomas, todos os genomas existentes

#Fiz isto originalmente...
#        for value in range(len(dic_trackback[key])): #ir às proteinas desses genomas
#            if dic_trackback[key][value] not in proteins: #se essas ainda não estiverem na lista de proteinas
#                proteins.append(dic_trackback[key][value]) #colocar essas proteinas na lista de proteinas
    
    #Aqui ja tenho as listas de proteinas e genomas direitinhas! 
    
    #Passo2_ Inicializar a matriz só com 0's : nr de linhas = nr de genomas, nr de colunas = nr de proteinas
    matrix=[] 
    for i in range(len(genomes)):
        matrix.append([0]) #numero de linhas igual ao numero de genomas
        for j in range(len(proteins)-1):
            matrix[i].append(0)#numero de colunas igual ao numero de proteinas
            
    #Passo3
    #1a parte:
    #Percorrer as chaves do dic_clusters, ver a que genoma pertence cada um dos seus values (proteinas)
    #inserir o 1 na matriz

#Originalmente andava a mexer no dic_blast... e vi as keys e depois os values e fui a ver que genoma pertenciam   
#    for key in dic_blast:
#        genome_name=''
#        linha = None #qual vai ser a linha
#        coluna = None
#        for gen in dic_trackback:
#            if key in dic_trackback[gen]:# se a proteina estiver nos valores desse genoma
#                genome_name=gen #guardar numa variavel o nome desse genoma
#                break
#         #procurar na lista de genomas o genoma anterior e saber a sua posiçao
#        linha= genomes.index(genome_name)
#        coluna= proteins.index(key)
#        matrix[linha][coluna]=1
 
#        
#    #2a parte: Ir aos valores de cada proteina e fazer o mesmo (a parte mais dificil)
#        
#        for value in range(len(dic_blast[key])):
#            geno_name=''
#            lin=None
#            col=None
#            for geno in dic_trackback:
#                if dic_blast[key][value] in dic_trackback[geno]: 
#                    geno_name= geno
#                    break
#            lin=genomes.index(geno_name)
#            col=proteins.index(dic_blast[key][value])
#            matrix[lin][col]=1
#            
#            matrix[linha][col]=1 #ir ao genoma da key e por um 1 na coluna do value
#            matrix[lin][coluna]=1# ir ao genoma do value e por um 1 na coluna da key
    
    
    #Passo3
    #Ir às chaves do dic_clusters, ver os values e de que genoma vêm
    #brotar de um 1
    
    for key in dic_clusters:
        coluna= proteins.index(key)
        for value in range(len(dic_clusters[key])):
            genome_name=''
            for gen in dic_trackback:
                if dic_clusters[key][value] in dic_trackback[gen]:# se a proteina estiver nos valores desse genoma
                    genome_name = gen #guardar numa variavel o nome desse genoma
                    break
#            procurar na lista de genomas o genoma anterior e saber a sua posiçao
        
            linha= genomes.index(genome_name)
            
            matrix[linha][coluna]=1
    
    
#Imprimir a matriz

    file = open("matrix.csv","w")
    file.write('sep=,\n')
    file.write('Phages,') #escreve phages no primeiro quadrado
    csv_writer = csv.writer(file, delimiter=',')
    csv_writer.writerow(proteins) #escrever so uma linha com o nome das proteinas
    
    
    for line in range(len(matrix)):
        file.write(genomes[line] + ',') # para escrever o nome do genoma
        csv_writer.writerow(matrix[line])
        #print(matrix[line])
    
    file.write('\n')
    for key in dic_clusters:
        if key.startswith('Cluster'):    
            file.write(key + '->' + '; '.join(dic_clusters[key])+'\n')
        
    file.close()
        

if __name__== '__main__':
    file=open('genomes.txt', 'w') #abrir um ficheiro para ser a database
    file.close()
    acNumbers= input('Write the acession numbers (separated by commas with no space):') #'MF033347,MF033348,MF033349,KY363465,NC_024124'
    dic_trackback= seq_NCBI(acNumbers) #receber os varios acNumbers, obter as suas proteinas, e colocar no ficheiro anterior
    clusters = blast()
    matrix(dic_trackback,  clusters)
   #MF033347,MF033348,MF033349,KY363465,KY471386.1,KY442063.1,MF428478.1,MF428479.1,MG656408.1,MF405094.1
