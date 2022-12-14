# Import parts of Biopython
from Bio import SeqIO
# import anytree
# import collections
import os
import numpy as np


# the function takes a fasta file and extracts the genome as a string
def extractFasta(path_to_file):
    sequence = ''
    with open(path_to_file, mode='r') as handle:
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            description = record.description
            sequence = record.seq

            # Example: adapt to extract features you are interested in
            print('----------------------------------------------------------')
            print('Processing the record {}:'.format(identifier))
            print('Its description is: \n{}'.format(description))
            amount_of_nucleotides = len(sequence)
            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides))

    sequence = str(sequence)
    print(type(sequence))
    # print(sequence)
    return sequence[0:49]


# the function cuts the sequence to kmers of length k and step
# returns a list with all the sequences
def kmers(k, step, sequence):
    # string_name[start:end:step]
    list1 = []
    for index in range(0, len(sequence), step):
        # print(sequence[index:index+k:step],'\n')
        list1.append(sequence[index:index + k:1])
    return list1

#
# !!!extend function to work on comparing the new genome to deeper lineages - done?
# (not use the basic genome, but the one that is the closest to him on the phylogenetic tree)
# the function compares the new genome to the base covid19 and saves the differences
def find_differntial(kmers_list, basic, families):
    min_dif = basic
    min_index = 0
    temp_arr = []
    for lineage in range(0, len(families)):  # iterating over he lineages (pages) - rows of the "families" array
        min_index = 0
        temp_arr = []
        for position in range(0, len(families[lineage])):  # iterating over the kmers in the lineage - cols of the "families" array
            if (families[lineage][position] != kmers_list[position]):
                #np.insert(families,position,families[lineage][position])
                temp_arr.append(kmers_list[position])
        #after comparing the given genome to the current lineage
        #checking if the current differntial is the smallest one
        #if so, saving the location in order to insert the genome there at the end
        print(temp_arr)
        if(temp_arr != [] and len(temp_arr)<=len(min_dif)):
            min_dif = temp_arr
            min_index = lineage
    if(min_dif != basic):
        families.insert(min_index+1,min_dif)

#exporting the database to a file in the following format:
#lineage1\n
#lineage2\n
#:
"""
def exportToFile(families):
    for lineage in range(0, len(families)):
        lines = families(lineage)
        with open('readme.txt', 'w') as f:
            for line in lines:
                f.write(line)
            f.write('\n')
    #f.close()

"""


# דאטה בייס ראשון
# הפונקציה מקבלת גנומים ממשפחות שונות ומוצאת את הקמרים המייצגים של כל משפחה ושומרת רק אותם במבנה הנתונים (את המקטעים שמשותפים למספר משפחות צריך לזרוק)
# את מספר המשפחות שמעליו נקבע שקיימר הוא גנרי וניתן לזרוק אותו נקבע כפרמטר ונשחק איתו עד שנמצא ערך טוב מספיק  -  ערך סף
# מטרת הדאטה בייס היא לתקן את בעיית היישור (אלינמנט) על ידי שנסווג גנום חדש למשפחה כלשהי על ידי בדיקה אם יש לו את הגנומים שמייצגים את המשפחה הזו
#  אם ניתקל במקרה שסיווגנו גנום  שמתאים ליותר ממשפחה אחת נבדוק האם הקיימרים שהגדרנו כמייצגים הם אכן כאלה או שהם גנרים וצריך לזרוק אותם ונבדוק גם האם בחרנו ערך סף טוב מספיק
# הפונקציה הזו תיצור את מבנה הנתונים הראשון , כלומר תשווה בין כל המשפחות ותמצא את הגנומים המייצגים של כל משפחה ותשמור אותם במבנה.
def database1_build(database1_built, rep_threshold): #rep_threshold - הסף שצריך לעבור בשביל לומר שהקיימר מייצג את המשפחה
    #database1 is the collection if sets we get from Build function below
    # what we want to receive in the function in order to work on in this function -
    # either raw fasta files of lineages and we will arrange them an find representatives -
    # receive from Build function the processed one
    #find representative kmers
    to_remove = set()
    count = 0
    for family in database1_built:
        for kmer in family:
            for index in range(0,len(database1_built)):
                if((database1_built[index]!=family) and (kmer in database1_built[index])):
                    count = count+1
            if(count>=rep_threshold):
                to_remove.add(kmer)
            count = 0
    #removing the generic kmers
    for removed in to_remove:
        for index2 in range(0, len(database1_built)):
            database1_built[index2].remove(removed)



#input: list of fasta files that represents a lineage
#output: array of sets of the lineage
#arranging the list of fasta files as kmers in sets
#a set contains the kmers of a lineage
def Build(list_of_fasta_files , k , step): #zuher's idea to implement!!!!!!!!!!!!!!!!!!!!!

    lineage = {} #empty set
    for i in list_of_fasta_files:
        lineage.add(set(kmers(k,step,i)))


# בהינתן מבנה הנתונים הראשון שקיים , הפונקציה תקבל גנום חדש ותקבע לאיזו משפחה הוא שייך (לפי היטים, כלומר פגיעות ומערך מונים)
# יש פה 2 מערכי מונים שונים !
# לגנום חדש שיש לסווג יש מערך מונים משלו, שמכיל כמה קיימרים מייצגים מכל משפחה יש בו
# אם הגנום החדש מכיל קמייר מייצג ששייך למשפחה מסוימת נגדיל את הקאונטר של הגנום ב1, כלומר נרצה למנות כמה קיימרים מייצגים של משפחה מסוימת קיימים בגנום שאנחנו רוצים לסווג
# נשים לב , על מנת לקבוע שהגנום אכן שייך למשפחה נגדיר ערך סף של מספר הקיימרים המייצגים שצריך שהוא יכיל על מנת שנוכל לאמר שהוא אכן במשפחה
# כדי לסווג גנום ששייך למשפחה , נגדיר ערך סף למספר הקיימרים המייצגים שהגנום צריך להכיל כדי להשתייך למשפחה מסוימת . אם הוא עובר את ערך הסף הזה - הוא שייך למשפחה וזו פגיעה
# את הפגיעות יש לשמור במערך מונים נוסף, כל תא במערך מייצג משפחה (המערך מאותחל לאפס) , פגיעה היא העלאת הערך בתא מסוים ב1 , בסוף התהליך עוברים על המערך , אם יש מספר תאים שבהם יש פגיעה , הגנום יסווג לתא עם הקאונטר הגבוה ביותר (שבו יש מספר גבוה של פגיעות)
# אם יש מספר תאים שבהם יש מספר גבוה של פגיעות , ככל הנראה מדובר בווריאנט חדש ששהוא המך השושלת של אחת מהמשפחות הללו

#איך מטפלים במקרה שבו עוברים את סף הסיווג ביותר בלפחות שני קאונטרים (כלומר יש סיווג ליותר ממשפחה אחת )
# rep_threshold המשמעות היא שהסף או שהוא נמוך מית=די או שהסף של הפונקציה הראשונה לא טוב , מכיוון שאם הסף של הפונקציה הראשונה נמוך מידי זה אומר שלקחנו קימר גנרי והגדרנו שהוא מייצג של יותר ממשפחה אחת לכן צריך לזרוק את הקיימר הזה , ולהעלות את הסף של הפונקציה הראשונה

def database1_classify(database1 , genome_to_classify, hit_threshold): #hit_threshold - הסף שצריך לעבור בשביל שגנום יסווג למשפחה כלשהי
    # we assume the database is an existing array
    #we got a genome that is already cut to kmers
    #מונה שבו כל תא מייצג משפחה וכאשר יש פגיעה בין הגנום לאחד המייצגים נעלה ב-1
    hit_counter = [0]*len(database1)
    #l=[0]*n
    # go through all the data base
    for kmer in genome_to_classify:
        #now we are in a specific family
        # אנחנו רוצים להשוות בין קיימר ב genome_to_classify לבין רשימת הקיימרים המייצגים שיש למשפחה
        # ההשוואה עצמה
        # genome_to_classify עושים השוואה בין הקיימר המייצג הנוכחי לבין כל הקיימרים שיש ל
        for family in range(0,len(database1)):
            if(kmer in database1[family]):
                #יש פגיעה נעלה את הקאונטר המתאים
                hit_counter[family] +=1  # זה לא ירוץ כיfamily  לא מספר אבל הכוונה זה לתא במערך המונים שמתאים למשפחה
    print(hit_counter)
        # עברנו על כל הקיימרים המייצגים שיש במשפחה
    # עברנו על כל המשפחות , נצפה לקבל מערך שיש בו תאים עם מספר גבוה/ נמוך של פגיעות
    # אם למשל במערך  hit_counter במקום ה0 (וזה תואם למשפחה A) יש ערך 3 , זה אומר שבגנום שאנחנו genome_to_classify יש 3 קיימרים שונים שתואמים ל3 קיימרים מייצגים במשפחה A
    # מבחינת סיווג , נחזיר את הaccession
    #כלומר התז של הlineage שהוא מסווג אליו
    
if __name__ == '__main__':
    # input_file =
    # output_file =
    # fasta_sequenses = SeqIO.parse(open(input_file),'fasta')
    # with open(output_file) as out_file:
    #     for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         new_sequence = some_function(sequence)
    #         write_fasta(out_file)

    # Install the biopython library (if not already installed)

    families = []
    path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/fasta/EPI_ISL_402124.fasta'  # <--- substitute by your local path
    sequence = extractFasta(path_to_file)
    kmers_list = kmers(4, 4, sequence)
    print(kmers_list)
    # print(type(kmers_list[7]))
    families.append(kmers_list)
    # families=np.array(families)

    directory = 'fasta'
    for filename in os.scandir(directory):
        if filename.is_file() and filename.path != 'fasta\EPI_ISL_402124.fasta':
            print(filename.path)
            # File path to your FASTA file
            path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/' + filename.path  # <--- substitute by your local path
            sequence = extractFasta(path_to_file)
            kmers_list = kmers(4, 4, sequence)
            # saving the differences between the new genome and the basic genome - !!compare new genome to the closest one in the lineage that is already in the database
            find_differntial(kmers_list, families[0], families)
            #print(kmers_list)
            # print(type(kmers_list[7]))
            #families.append(kmers_list)

    #print(families[0])  # the first genome is the reference basic covid19
    print(len(families))  # num of rows
    #print(len(families[0]))  # num of columns
    print('-------------------------------------------')
    print('find representative kmers')
    families2 = [{'TTTC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATTA', 'CAGG', 'AAAC', 'TTCC', 'TAAC'},{'TTTC', 'CAAC', 'TAAA', 'CTTT', 'TTGT', 'AAAA', 'AGAT', 'CTGT', 'TCTG', 'TCGA', 'ACTT'},{'TTTC', 'CAAC', 'CGAT', 'GAAC', 'CTCT', 'TTTA', 'AAAT', 'TGTA', 'CTTT', 'AAAC', 'CCAA', 'C'}]
    print(families2)
#   for index in range(0,len(families)):
#       families[index] = set(families[index])
#   print(families)
#   print(type(families[0]))
    print('********')
    database1_build(families2,2)
    print(families2)
    #print("Try writing to file")
    #exportToFile(families)
    print('-------------------------------------------')
    print("classify the genom to a specific family")
    database1_classify(families2,{'TATC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATGA', 'CAGG', 'ACAC', 'TTCC', 'TAAC'},0)
    d

