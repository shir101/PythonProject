# Import parts of Biopython
from Bio import SeqIO
import os
import numpy as np
from DataCollector import Covid19DataPortalAccessionFetcher as cdc


# input: list of fasta files that represents a lineage
# output: array of sets of the lineage
# arranging the list of fasta files as kmers in sets
# a set contains the kmers of a lineage
# will convert list of fasta files to list of kmers
def Build(list_of_fasta_files_accessions, k, step):  # zuher's idea to implement!!!!!!!!!!!!!!!!!!!!!
    lineage = []  # empty array of sets
    directory = 'raw'
    for filename in os.scandir(directory):
        if filename.is_file():
            acc_set = {}
            #print(filename.path)
            # File path to your FASTA file
            path_to_file = 'Home/Desktop/ViralDataCollector-main/data/raw' + filename.path  # <--- substitute by your local path
            sequence = extractFasta(path_to_file)
            acc_set = set(kmers(k, step, sequence))
            # kmers_list = kmers(k, step, sequence)
            lineage.add(acc_set)

    # TODO: takes the fasta files and converts them to an array that is usable for the first database - check this function!!!!!!
    # for i in list_of_fasta_files:
    #    lineage.add(set(kmers(k,step,i)))


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
        for position in range(0, len(
                families[lineage])):  # iterating over the kmers in the lineage - cols of the "families" array
            if (families[lineage][position] != kmers_list[position]):
                # np.insert(families,position,families[lineage][position])
                temp_arr.append(kmers_list[position])
        # after comparing the given genome to the current lineage
        # checking if the current differntial is the smallest one
        # if so, saving the location in order to insert the genome there at the end
        #print(temp_arr)
        if (temp_arr != [] and len(temp_arr) <= len(min_dif)):
            min_dif = temp_arr
            min_index = lineage
    if (min_dif != basic):
        families.insert(min_index + 1, min_dif)


# ----------------------------------------------------------------------------------


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
            #print('----------------------------------------------------------')
            #print('Processing the record {}:'.format(identifier))
            #print('Its description is: \n{}'.format(description))
            amount_of_nucleotides = len(sequence)
            #print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides))

    sequence = str(sequence)
    #print(type(sequence))
    # print(sequence)
    return sequence[0:200]


# the function cuts the sequence to kmers of length k and step
# returns a list with all the sequences
def kmers(k, step, sequence):
    # string_name[start:end:step]
    list1 = []
    for index in range(0, len(sequence), step):
        # print(sequence[index:index+k:step],'\n')
        # this if is for removing the extra kmers that are created when reaching the end of the sequence
        # when the length is smaller than k it means we reached the end, and the small sequences created
        # are already covered in the last kmer
        if len(sequence[index:index + k:1])<k:
            pass
        else:
            list1.append(sequence[index:index + k:1])
    return set(list1)


# הפונקציה מקבלת את קבצי ה-accessions של lineage אחד
# מוצאת את הקיימרים המייצגים של ה-lineage מתוך ה-accessions שקיבלנו
def lineage_reps_from_its_accessions(dict_of_acc, rep_threshold):
    # acc -
    # rep_threshold - if above the threshold than the kmer is representative
    # returns a set of the representative kmers of the lineage from the accessions
    to_keep = set()
    count = 1
    for accession in dict_of_acc:
        for kmer in dict_of_acc[accession]:
            search_value = kmer
            matches = {k: v for k, v in dict_of_acc.items() if search_value in v}
            # print(matches)
            count = len(matches)
            # for check in dict_of_acc:
            #   if ((check != accession) and (kmer in check)):
            # count = count + 1
            if (count > rep_threshold):  # if lower than threshold than not representing so we throw it
                to_keep.add(kmer)
            count = 1
    return to_keep


# דאטה בייס ראשון
# מסתכלים על lineaage יחיד כלשהו, ולוקחים הרבה קבצי פאסטה שמייצגים accessions שלו
# נעבור על ה-accessions השונים ונמצא את הקיימרים המייצגים של כל lineage
# נעשה זאת עבור כל lineage בנפרד (לא לערבב בטעות בין lineages שונים - זה של דאטה בייס 2)
# הערה חשובה: לlineage אחד יש כמה accessions ! , למשל ל B.1.1.7 יש OE996776 , MZ166750 וכו
def database1_build(dict_of_acc_of_a_lineage):  # fit to Zuher's interface
    # database will be a single lineage and its accessions
    lineage_name = dict_of_acc_of_a_lineage
    # עוברים על כל ה lineages, ועבור כל אחד מהם, מוציאם את המייצגים מתוך ה-accessions
    # acc = [] #array of sets (like families)
    rep_threshold = 0  # need to either get it or findn it using optimizations!!!!
    # for ...
    #    set_of_reps_for_lineage = lineage_reps_from_its_accessions(acc, rep_threshold)
    #    database[lineage] = set_of_reps_for_lineage
    # TODO: insert the set of representatives of each lineage to the database

    # finish comparing one to all the other

    # returning the representative k-mers of a single lineage
    return lineage_reps_from_its_accessions(dict_of_acc_of_a_lineage, rep_threshold)

    # use this function to optimize the threshold
    # database1 comprehensionbe done in the main


# דאטה בייס שני
# הפונקציה מקבלת גנומים ממשפחות שונות ומוצאת את הקמרים המייצגים של כל משפחה ושומרת רק אותם במבנה הנתונים (את המקטעים שמשותפים למספר משפחות צריך לזרוק)
# את מספר המשפחות שמעליו נקבע שקיימר הוא גנרי וניתן לזרוק אותו נקבע כפרמטר ונשחק איתו עד שנמצא ערך טוב מספיק  -  ערך סף
# מטרת הדאטה בייס היא לתקן את בעיית היישור (אלינמנט) על ידי שנסווג גנום חדש למשפחה כלשהי על ידי בדיקה אם יש לו את הגנומים שמייצגים את המשפחה הזו
#  אם ניתקל במקרה שסיווגנו גנום  שמתאים ליותר ממשפחה אחת נבדוק האם הקיימרים שהגדרנו כמייצגים הם אכן כאלה או שהם גנרים וצריך לזרוק אותם ונבדוק גם האם בחרנו ערך סף טוב מספיק
# הפונקציה הזו תיצור את מבנה הנתונים השני , כלומר תשווה בין כל המשפחות ותמצא את הגנומים המייצגים של כל משפחה ותשמור אותם במבנה.
# השוואה בין קיימרים מייצגים של lineages וזריקת המיותרים
def database2_build(database2_built,
                    dif_threshold):  # dif_threshold - הסף שצריך לעבור בשביל לומר שהקיימר מייצג את המשפחה, מה שעובר אותו אז צריך להסיר כי גנרי מדי
    # database2 is the collection of sets we get from Build function below
    # what we want to receive in the function in order to work on in this function -
    # either raw fasta files of lineages and we will arrange them an find representatives -
    # receive from Build function the processed one
    # find representative kmers
    to_remove = set()
    count = 0
    for family in database2_built:
        for kmer in database2_built[family]:
            for index in database2_built:
                if ((database2_built[index] != database2_built[family]) and (kmer in database2_built[index])):
                    count = count + 1
            if (count >= dif_threshold):  # if bigger than threshold than too generic and not representing
                to_remove.add(kmer)
            count = 0
    # removing the generic kmers
    for index2 in database2_built:
        for removed in to_remove:
            temp = database2_built[index2]
            temp.discard(removed)
            database2_built[index2] = temp
    return database2_built



# בהינתן מבנה הנתונים השני שקיים , הפונקציה תקבל גנום חדש ותקבע לאיזו משפחה הוא שייך (לפי היטים, כלומר פגיעות ומערך מונים)
# יש פה 2 מערכי מונים שונים !
# לגנום חדש שיש לסווג יש מערך מונים משלו, שמכיל כמה קיימרים מייצגים מכל משפחה יש בו
# אם הגנום החדש מכיל קמייר מייצג ששייך למשפחה מסוימת נגדיל את הקאונטר של הגנום ב1, כלומר נרצה למנות כמה קיימרים מייצגים של משפחה מסוימת קיימים בגנום שאנחנו רוצים לסווג
# נשים לב , על מנת לקבוע שהגנום אכן שייך למשפחה נגדיר ערך סף של מספר הקיימרים המייצגים שצריך שהוא יכיל על מנת שנוכל לאמר שהוא אכן במשפחה
# כדי לסווג גנום ששייך למשפחה , נגדיר ערך סף למספר הקיימרים המייצגים שהגנום צריך להכיל כדי להשתייך למשפחה מסוימת . אם הוא עובר את ערך הסף הזה - הוא שייך למשפחה וזו פגיעה
# את הפגיעות יש לשמור במערך מונים נוסף, כל תא במערך מייצג משפחה (המערך מאותחל לאפס) , פגיעה היא העלאת הערך בתא מסוים ב1 , בסוף התהליך עוברים על המערך , אם יש מספר תאים שבהם יש פגיעה , הגנום יסווג לתא עם הקאונטר הגבוה ביותר (שבו יש מספר גבוה של פגיעות)
# אם יש מספר תאים שבהם יש מספר גבוה של פגיעות , ככל הנראה מדובר בווריאנט חדש ששהוא המך השושלת של אחת מהמשפחות הללו

# איך מטפלים במקרה שבו עוברים את סף הסיווג ביותר בלפחות שני קאונטרים (כלומר יש סיווג ליותר ממשפחה אחת )
# rep_threshold המשמעות היא שהסף או שהוא נמוך מידי או שהסף של הפונקציה הראשונה לא טוב , מכיוון שאם הסף של הפונקציה הראשונה נמוך מידי זה אומר שלקחנו קימר גנרי והגדרנו שהוא מייצג של יותר ממשפחה אחת לכן צריך לזרוק את הקיימר הזה , ולהעלות את הסף של הפונקציה הראשונה

def database2_classify(database2, genome_to_classify,
                       hit_threshold):  # hit_threshold - הסף שצריך לעבור בשביל שגנום יסווג למשפחה כלשהי
    # we assume the database is an existing array
    # we got a genome that is already cut to kmers
    # מונה שבו כל תא מייצג משפחה וכאשר יש פגיעה בין הגנום לאחד המייצגים נעלה ב-1
    hit_counter = [0] * len(database2)
    # l=[0]*n
    # go through all the data base
    for kmer in genome_to_classify:
        # now we are in a specific family
        # אנחנו רוצים להשוות בין קיימר ב genome_to_classify לבין רשימת הקיימרים המייצגים שיש למשפחה
        # ההשוואה עצמה
        # genome_to_classify עושים השוואה בין הקיימר המייצג הנוכחי לבין כל הקיימרים שיש ל
        for family in range(0, len(database2)):
            if (kmer in database2[family]):
                # יש פגיעה נעלה את הקאונטר המתאים
                hit_counter[family] += 1  # זה לא ירוץ כיfamily  לא מספר אבל הכוונה זה לתא במערך המונים שמתאים למשפחה
    print(hit_counter)
    lineage_of_classification = classify_hits_by_threshold(database2, hit_counter,
                                                           hit_threshold)  # TODO: return the lineage by using a dictionary

    return lineage_of_classification
    # עברנו על כל הקיימרים המייצגים שיש במשפחה
    # עברנו על כל המשפחות , נצפה לקבל מערך שיש בו תאים עם מספר גבוה/ נמוך של פגיעות
    # אם למשל במערך  hit_counter במקום ה0 (וזה תואם למשפחה A) יש ערך 3 , זה אומר שבגנום שאנחנו genome_to_classify יש 3 קיימרים שונים שתואמים ל3 קיימרים מייצגים במשפחה A
    # מבחינת סיווג , נחזיר את הaccession
    # כלומר התז של הlineage שהוא מסווג אליו


# מתודת עזר למתודת ה-classify
# מבצעת את הסיווג המתאים לפי ערך הסף ומתחשבת במקרים שונים
# תבצע את החזרת ה-accession המתאים לפי הסיווג והסף
def classify_hits_by_threshold(database1, hit_counter, hit_threshold):
    max_hits = max(hit_counter)  # max value of hits
    get_indexes = [i for i, x in enumerate(hit_counter) if
                   x == max_hits]  # all the indexes that contain the max number of hits
    # TODO: finish this method


if __name__ == '__main__':
    #    families = []
    #   path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/fasta/EPI_ISL_402124.fasta'  # <--- substitute by your local path
    #    sequence = extractFasta(path_to_file)
    #    kmers_list = kmers(4, 4, sequence)
    #   print(kmers_list)
    # print(type(kmers_list[7]))
    #    families.append(kmers_list)
    # families=np.array(families)

    #    directory = 'fasta'
    #    for filename in os.scandir(directory):
    #        if filename.is_file() and filename.path != 'fasta\EPI_ISL_402124.fasta':
    #            print(filename.path)
    # File path to your FASTA file
    #            path_to_file = 'C:/Users/User/PycharmProjects/DatabaseBuild/' + filename.path  # <--- substitute by your local path
    #            sequence = extractFasta(path_to_file)
    #            kmers_list = kmers(4, 4, sequence)
    # saving the differences between the new genome and the basic genome - !!compare new genome to the closest one in the lineage that is already in the database
    #            find_differntial(kmers_list, families[0], families)
    # print(kmers_list)
    # print(type(kmers_list[7]))
    # families.append(kmers_list)

    # print(families[0])  # the first genome is the reference basic covid19
    #    print(len(families))  # num of rows
    # print(len(families[0]))  # num of columns
    #    print('-------------------------------------------')
    #    print('find representative kmers')
    # families2 = [{'TTTC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATTA', 'CAGG', 'AAAC', 'TTCC', 'TAAC'},{'TTTC', 'CAAC', 'TAAA', 'CTTT', 'TTGT', 'AAAA', 'AGAT', 'CTGT', 'TCTG', 'TCGA', 'ACTT'},{'TTTC', 'CAAC', 'CGAT', 'GAAC', 'CTCT', 'TTTA', 'AAAT', 'TGTA', 'CTTT', 'AAAC', 'CCAA', 'C'}]
    # print(families2)
    #   for index in range(0,len(families)):
    #       families[index] = set(families[index])
    #   print(families)
    #   print(type(families[0]))
    # print('********')
    # database2_build(families2,2)
    # print(families2)
    # print("Try writing to file")
    # exportToFile(families)
    # print('-------------------------------------------')
    # print("classify the genome to a specific family")
    # database2_classify(families2,{'TATC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATGA', 'CAGG', 'ACAC', 'TTCC', 'TAAC'},0)
    # print("representing kmers of lineage by accessions")
    # lineage_reps = lineage_reps_from_its_accessions(families2, 3)
    # print(lineage_reps)
    # classify_hits_by_threshold(families2, [8,0,1,8,3,8,8], 0)

    ####-------------------------------ignore-above------------------------------------

    # This part of code is responsible for the extraction of fasta files from the data/raw directory!!!!!!! - important
    # print(os.getcwd())
    # print("files")
    # directory  = './data/raw'  # <--- substitute by your local path
    # for filename in os.listdir(directory):
    #    f = os.path.join(directory, filename)
    # checking if it is a file
    #    if os.path.isfile(f):
    #        print(f)

    #dc = cdc() # instance of the data collector
    #print("dc cdc successful!!!!!")
    #avail_lineages = dc.getLineagesList() ####
    #avail_lineages = ['B.1.1', 'B.1.1.285', 'B.12', 'B.1.1.284', 'B.1.1.214']
    # print(avail_lineages)
    
    #building database1 for all the lineages 
    database1 = {}  # global dictionary for database 1
    #for lineage in avail_lineages[1:5]:
        #accessions = dc.getAccessionsByLineage(lineage) #### take accessions of a certain lineage
        #dc.downloadAccessions(accessions)
        # accessions = []
        # print(accessions)

    directory = './data/'  # <--- substitute by your local path

    for dirname in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, dirname)):
            print(dirname)

        try:
            if dirname=='.ipynb_checkpoints' or dirname=='raw':
                continue
            else:
                dict = {}
                for filename in os.listdir(os.path.join(directory, dirname)):
                    f = os.path.join(os.path.join(directory, dirname), filename)
                    # checking if it is a file
                    if os.path.isfile(f):
                        sequence = extractFasta(f)
                        kmers_list = kmers(32, 1, sequence)
                        #print(kmers_list)
                        dict[filename] = kmers_list
                #print(dict)
                #print(database1_build(dict))

                reps_of_lineage = lineage_reps_from_its_accessions(dict, 13)
                print("**************************")
                #print("We can see here a dictionary of 1 lineage and its representative genomes")
                database1.update({dirname: reps_of_lineage})
                #print(dict)

        except FileNotFoundError:
            pass

        #delete the files in the data/raw of the lineage after inserting it to database1
        #for filename in os.listdir(directory):
            #f = os.path.join(directory, filename)
            # checking if it is a file
            #if os.path.isfile(f):
                #os.remove(f)

    print(database1)
    print('------------------------database-2 -------------------------------------')
    print(database2_build(database1,2))
    #print(database1)
