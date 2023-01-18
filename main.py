# Import parts of Biopython
from Bio import SeqIO
import os
import numpy as np
from DataCollector import Covid19DataPortalAccessionFetcher as cdc




#input: list of fasta files that represents a lineage
#output: array of sets of the lineage
#arranging the list of fasta files as kmers in sets
#a set contains the kmers of a lineage
#will convert list of fasta files to list of kmers
def Build(list_of_fasta_files_accessions , k , step): #zuher's idea to implement!!!!!!!!!!!!!!!!!!!!!
    lineage = []  # empty array of sets
    directory = 'raw'
    for filename in os.scandir(directory):
        if filename.is_file():
            acc_set = {}
            print(filename.path)
            # File path to your FASTA file
            path_to_file = 'Home/Desktop/ViralDataCollector-main/data/raw' + filename.path  # <--- substitute by your local path
            sequence = extractFasta(path_to_file)
            acc_set = set(kmers(k, step, sequence))
            #kmers_list = kmers(k, step, sequence)
            lineage.add(acc_set)

    #TODO: takes the fasta files and converts them to an array that is usable for the first database - check this function!!!!!!
    #for i in list_of_fasta_files:
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

print(

#----------------------------------------------------------------------------------


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
    return set(list1)






#הפונקציה מקבלת את קבצי ה-accessions של lineage אחד
#מוצאת את הקיימרים המייצגים של ה-lineage מתוך ה-accessions שקיבלנו
def lineage_reps_from_its_accessions(dict_of_acc, rep_threshold):
    # acc -
    # rep_threshold - if above the threshold than the kmer is representative
    # returns a set of the representative kmers of the lineage from the accessions
    to_keep = set()
    count = 1
    for accession in dict_of_acc:
        for kmer in dict_of_acc[accession]:
            search_value =  kmer
            matches = {k: v for k, v in dict_of_acc.items() if search_value in v}
            #print(matches)
            count = len(matches)
            #for check in dict_of_acc:
             #   if ((check != accession) and (kmer in check)):
                   # count = count + 1
            if (count > rep_threshold):  # if lower than threshold than not representing so we throw it
                to_keep.add(kmer)
            count = 1
    return to_keep


#דאטה בייס ראשון
#מסתכלים על lineaage יחיד כלשהו, ולוקחים הרבה קבצי פאסטה שמייצגים accessions שלו
#נעבור על ה-accessions השונים ונמצא את הקיימרים המייצגים של כל lineage
#נעשה זאת עבור כל lineage בנפרד (לא לערבב בטעות בין lineages שונים - זה של דאטה בייס 2)
#הערה חשובה: לlineage אחד יש כמה accessions ! , למשל ל B.1.1.7 יש OE996776 , MZ166750 וכו
def database1_build(dict_of_acc_of_a_lineage):  #fit to Zuher's interface
    # database will be a single lineage and its accessions
    lineage_name = dict_of_acc_of_a_lineage
    # עוברים על כל ה lineages, ועבור כל אחד מהם, מוציאם את המייצגים מתוך ה-accessions
    #acc = [] #array of sets (like families)
    rep_threshold = 0 #need to either get it or findn it using optimizations!!!!
    #for ...
    #    set_of_reps_for_lineage = lineage_reps_from_its_accessions(acc, rep_threshold)
    #    database[lineage] = set_of_reps_for_lineage
    #TODO: insert the set of representatives of each lineage to the database
    
    
    # finish comparing one to all the other
        
    #returning the representative k-mers of a single lineage
    return lineage_reps_from_its_accessions(dict_of_acc_of_a_lineage, rep_threshold)
                
    #use this function to optimize the threshold
    #database1 comprehensionbe done in the main
    



# דאטה בייס שני
# הפונקציה מקבלת גנומים ממשפחות שונות ומוצאת את הקמרים המייצגים של כל משפחה ושומרת רק אותם במבנה הנתונים (את המקטעים שמשותפים למספר משפחות צריך לזרוק)
# את מספר המשפחות שמעליו נקבע שקיימר הוא גנרי וניתן לזרוק אותו נקבע כפרמטר ונשחק איתו עד שנמצא ערך טוב מספיק  -  ערך סף
# מטרת הדאטה בייס היא לתקן את בעיית היישור (אלינמנט) על ידי שנסווג גנום חדש למשפחה כלשהי על ידי בדיקה אם יש לו את הגנומים שמייצגים את המשפחה הזו
#  אם ניתקל במקרה שסיווגנו גנום  שמתאים ליותר ממשפחה אחת נבדוק האם הקיימרים שהגדרנו כמייצגים הם אכן כאלה או שהם גנרים וצריך לזרוק אותם ונבדוק גם האם בחרנו ערך סף טוב מספיק
# הפונקציה הזו תיצור את מבנה הנתונים השני , כלומר תשווה בין כל המשפחות ותמצא את הגנומים המייצגים של כל משפחה ותשמור אותם במבנה.
#השוואה בין קיימרים מייצגים של lineages וזריקת המיותרים
def database2_build(database2_built, dif_threshold): #dif_threshold - הסף שצריך לעבור בשביל לומר שהקיימר מייצג את המשפחה, מה שעובר אותו אז צריך להסיר כי גנרי מדי
    #database2 is the collection of sets we get from Build function below
    # what we want to receive in the function in order to work on in this function -
    # either raw fasta files of lineages and we will arrange them an find representatives -
    # receive from Build function the processed one
    #find representative kmers
    to_remove = set()
    count = 0
    for family in database2_built:
        for kmer in family:
            for index in range(0,len(database2_built)):
                if((database2_built[index]!=family) and (kmer in database2_built[index])):
                    count = count+1
            if(count>=dif_threshold): #if bigger than threshold than too generic and not representing
                to_remove.add(kmer)
            count = 0
    #removing the generic kmers
    for removed in to_remove:
        for index2 in range(0, len(database2_built)):
            database2_built[index2].remove(removed)



# בהינתן מבנה הנתונים השני שקיים , הפונקציה תקבל גנום חדש ותקבע לאיזו משפחה הוא שייך (לפי היטים, כלומר פגיעות ומערך מונים)
# יש פה 2 מערכי מונים שונים !
# לגנום חדש שיש לסווג יש מערך מונים משלו, שמכיל כמה קיימרים מייצגים מכל משפחה יש בו
# אם הגנום החדש מכיל קמייר מייצג ששייך למשפחה מסוימת נגדיל את הקאונטר של הגנום ב1, כלומר נרצה למנות כמה קיימרים מייצגים של משפחה מסוימת קיימים בגנום שאנחנו רוצים לסווג
# נשים לב , על מנת לקבוע שהגנום אכן שייך למשפחה נגדיר ערך סף של מספר הקיימרים המייצגים שצריך שהוא יכיל על מנת שנוכל לאמר שהוא אכן במשפחה
# כדי לסווג גנום ששייך למשפחה , נגדיר ערך סף למספר הקיימרים המייצגים שהגנום צריך להכיל כדי להשתייך למשפחה מסוימת . אם הוא עובר את ערך הסף הזה - הוא שייך למשפחה וזו פגיעה
# את הפגיעות יש לשמור במערך מונים נוסף, כל תא במערך מייצג משפחה (המערך מאותחל לאפס) , פגיעה היא העלאת הערך בתא מסוים ב1 , בסוף התהליך עוברים על המערך , אם יש מספר תאים שבהם יש פגיעה , הגנום יסווג לתא עם הקאונטר הגבוה ביותר (שבו יש מספר גבוה של פגיעות)
# אם יש מספר תאים שבהם יש מספר גבוה של פגיעות , ככל הנראה מדובר בווריאנט חדש ששהוא המך השושלת של אחת מהמשפחות הללו

#איך מטפלים במקרה שבו עוברים את סף הסיווג ביותר בלפחות שני קאונטרים (כלומר יש סיווג ליותר ממשפחה אחת )
# rep_threshold המשמעות היא שהסף או שהוא נמוך מידי או שהסף של הפונקציה הראשונה לא טוב , מכיוון שאם הסף של הפונקציה הראשונה נמוך מידי זה אומר שלקחנו קימר גנרי והגדרנו שהוא מייצג של יותר ממשפחה אחת לכן צריך לזרוק את הקיימר הזה , ולהעלות את הסף של הפונקציה הראשונה

def database2_classify(database2 , genome_to_classify, hit_threshold): #hit_threshold - הסף שצריך לעבור בשביל שגנום יסווג למשפחה כלשהי
    # we assume the database is an existing array
    #we got a genome that is already cut to kmers
    #מונה שבו כל תא מייצג משפחה וכאשר יש פגיעה בין הגנום לאחד המייצגים נעלה ב-1
    hit_counter = [0]*len(database2)
    #l=[0]*n
    # go through all the data base
    for kmer in genome_to_classify:
        #now we are in a specific family
        # אנחנו רוצים להשוות בין קיימר ב genome_to_classify לבין רשימת הקיימרים המייצגים שיש למשפחה
        # ההשוואה עצמה
        # genome_to_classify עושים השוואה בין הקיימר המייצג הנוכחי לבין כל הקיימרים שיש ל
        for family in range(0,len(database2)):
            if(kmer in database2[family]):
                #יש פגיעה נעלה את הקאונטר המתאים
                hit_counter[family] +=1  # זה לא ירוץ כיfamily  לא מספר אבל הכוונה זה לתא במערך המונים שמתאים למשפחה
    print(hit_counter)
    lineage_of_classification = classify_hits_by_threshold(database2,hit_counter, hit_threshold) #TODO: return the lineage by using a dictionary

    return lineage_of_classification
    # עברנו על כל הקיימרים המייצגים שיש במשפחה
    # עברנו על כל המשפחות , נצפה לקבל מערך שיש בו תאים עם מספר גבוה/ נמוך של פגיעות
    # אם למשל במערך  hit_counter במקום ה0 (וזה תואם למשפחה A) יש ערך 3 , זה אומר שבגנום שאנחנו genome_to_classify יש 3 קיימרים שונים שתואמים ל3 קיימרים מייצגים במשפחה A
    # מבחינת סיווג , נחזיר את הaccession
    #כלומר התז של הlineage שהוא מסווג אליו

#מתודת עזר למתודת ה-classify
#מבצעת את הסיווג המתאים לפי ערך הסף ומתחשבת במקרים שונים
#תבצע את החזרת ה-accession המתאים לפי הסיווג והסף
def classify_hits_by_threshold(database1,hit_counter,hit_threshold):
    max_hits = max(hit_counter) #max value of hits
    get_indexes = [i for i, x in enumerate(hit_counter) if x==max_hits] #all the indexes that contain the max number of hits
    #TODO: finish this method




if __name__ == '__main__':
    # input_file =
    # output_file =
    # fasta_sequenses = SeqIO.parse(open(input_file),'fasta')
    # with open(output_file) as out_file:
    #     for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         new_sequence = some_function(sequence)
    #         write_fasta(out_file)


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
            #print(kmers_list)
            # print(type(kmers_list[7]))
            #families.append(kmers_list)

    #print(families[0])  # the first genome is the reference basic covid19
#    print(len(families))  # num of rows
    #print(len(families[0]))  # num of columns
#    print('-------------------------------------------')
#    print('find representative kmers')
    #families2 = [{'TTTC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATTA', 'CAGG', 'AAAC', 'TTCC', 'TAAC'},{'TTTC', 'CAAC', 'TAAA', 'CTTT', 'TTGT', 'AAAA', 'AGAT', 'CTGT', 'TCTG', 'TCGA', 'ACTT'},{'TTTC', 'CAAC', 'CGAT', 'GAAC', 'CTCT', 'TTTA', 'AAAT', 'TGTA', 'CTTT', 'AAAC', 'CCAA', 'C'}]
    #print(families2)
#   for index in range(0,len(families)):
#       families[index] = set(families[index])
#   print(families)
#   print(type(families[0]))
    #print('********')
    #database2_build(families2,2)
    #print(families2)
    #print("Try writing to file")
    #exportToFile(families)
    #print('-------------------------------------------')
    #print("classify the genome to a specific family")
    #database2_classify(families2,{'TATC', 'CAAC', 'GATC', 'AAGG', 'TACC', 'TTTA', 'T', 'ATGA', 'CAGG', 'ACAC', 'TTCC', 'TAAC'},0)
    #print("representing kmers of lineage by accessions")
    #lineage_reps = lineage_reps_from_its_accessions(families2, 3)
    #print(lineage_reps)
    #classify_hits_by_threshold(families2, [8,0,1,8,3,8,8], 0)
    
####-------------------------------ignore-above------------------------------------  

    # This part of code is responsible for the extraction of fasta files from the data/raw directory!!!!!!! - important
    #print(os.getcwd())
    #print("files")
    #directory  = './data/raw'  # <--- substitute by your local path
    #for filename in os.listdir(directory):
    #    f = os.path.join(directory, filename)
        # checking if it is a file
    #    if os.path.isfile(f):
    #        print(f)
            
    #dc = cdc() # instance of the data collector
    #avail_lineages = dc.getLineagesList() ####
    avail_lineages = ['B.1.1', 'B.1.1.285', 'B.12', 'B.1.1.284', 'B.1.1.214', 'B.1.346', 'B.1.1.45', 'B.1', 'B.1.2', 'B.1.1.485', 'B.1.1.174', 'B.1.1.7', '', 'R.1', 'B.1.1.101', 'B.1.1.17', 'B.1.1.482', 'B.1.1.372', 'AY.23', 'AY.29', 'BA.1.1.2', 'BA.1.1', 'AY.75.3', 'AY.85', 'BA.1.15', 'AY.29.1', 'AY.34.1', 'B.1.617.2', 'AY.29.2', 'AY.14', 'BA.2.3', 'BA.1.1.1', 'BA.1.17.2', 'BA.1.20', 'BA.1', 'BA.1.1.9', 'BA.2', 'BA.2.10', 'BA.2.29', 'BA.2.3.1', 'BA.2.24', 'BA.2.3.2', 'BA.2.10.1', 'B.1.177', 'B.1.258', 'B.1.1.213', 'B.1.13', 'B.1.391', 'B.1.104', 'B.1.439', 'B.1.93', 'B.1.177.82', 'B.1.1.164', 'B.1.1.217', 'B.1.1.296', 'AD.2', 'B.1.1.1', 'B.1.1.369', 'B.1.1.178', 'B.1.35', 'B.35', 'B.1.153', 'B.1.1.13', 'B.1.1.3', 'B.1.97', 'B.1.1.379', 'B.39', 'B.40', 'B.1.351', 'B.1.221', 'B.1.1.10', 'B.1.160.29', 'B.1.1.519', 'B.1.160', 'P.1', 'B.1.177.44', 'C.16', 'Q.4', 'B.1.258.17', 'B.1.160.16', 'B.1.160.30', 'B.1.177.62', 'B.1.36.1', 'Q.1', 'B.1.177.23', 'B.1.1.39', 'B.1.160.31', 'B', 'B.3', 'B.1.160.15', 'B.1.36.35', 'B.1.177.43', 'B.1.1.189', 'P.2', 'L.3', 'B.1.111', 'B.1.160.9', 'B.1.160.14', 'W.4', 'B.1.1.70', 'B.1.177.81', 'B.1.416.1', 'B.1.160.22', 'B.1.1.445', 'B.1.258.16', 'B.1.177.71', 'B.1.160.10', 'B.1.1.37', 'A.27', 'B.1.229', 'B.1.177.86', 'B.1.177.85', 'B.1.160.26', 'B.1.466', 'B.1.1.384', 'B.1.1.317', 'B.1.177.51', 'B.1.177.32', 'B.1.177.72', 'A.23.1', 'B.1.177.60', 'C.35', 'B.1.177.52', 'B.1.177.4', 'B.1.1.420', 'B.1.575', 'A.28', 'B.1.146', 'B.1.243', 'B.1.36', 'C.36', 'B.1.236', 'Y.1', 'B.1.588', 'B.1.214.2', 'B.1.36.31', 'B.1.480', 'B.1.1.219', 'B.1.427', 'B.1.177.77', 'B.1.1.222', 'B.1.214.3', 'B.1.142', 'B.1.1.89', 'B.1.201', 'B.1.465', 'B.1.76', 'B.1.456', 'B.1.1.274', 'B.1.1.8', 'B.1.356', 'A.2', 'B.1.36.29', 'B.1.470', 'B.1.619', 'B.1.466.1', 'A.21', 'B.1.177.28', 'B.1.131', 'B.1.258.14', 'B.1.525', 'B.1.177.15', 'B.1.1.232', 'B.1.367', 'B.1.147', 'B.1.474', 'B.1.128', 'B.1.332', 'B.1.12', 'B.1.406', 'B.1.160.32', 'B.1.149', 'B.1.160.12', 'B.52', 'B.1.105', 'B.1.1.194', 'B.31', 'B.58', 'B.1.1.198', 'B.1.1.197', 'B.27', 'B.1.379', 'B.1.157', 'B.1.1.71', 'B.23', 'B.3.1', 'B.29', 'B.1.167', 'B.1.383', 'B.1.1.171', 'B.1.1.166', 'B.1.1.220', 'B.28', 'B.1.39', 'B.1.1.419', 'B.1.81', 'B.1.1.270', 'B.1.117', 'B.1.198', 'B.1.1.378', 'B.1.1.279', 'B.4.7', 'B.1.250', 'A.2.2', 'B.1.1.277', 'B.1.254', 'B.1.1.291', 'B.1.1.75', 'B.1.1.275', 'B.1.1.391', 'B.1.1.257', 'B.10', 'B.1.151', 'B.1.1.83', 'B.1.1.43', 'B.1.1.254', 'B.1.22', 'B.1.1.4', 'B.33', 'B.11', 'B.1.223', 'B.61', 'B.1.9', 'A.1', 'B.1.1.304', 'B.1.1.41', 'B.1.1.5', 'B.1.1.208', 'B.55', 'B.1.91', 'B.1.1.25', 'B.1.321', 'B.1.165', 'B.1.1.323', 'B.57', 'B.1.491', 'B.1.1.168', 'B.45', 'B.34', 'B.32', 'B.26', 'A.5', 'A.2.3', 'B.1.222', 'B.1.559', 'B.1.256', 'B.1.248', 'B.1.215', 'B.1.220', 'B.1.557', 'B.1.1.119', 'A', 'B.1.1.395', 'C.27', 'B.1.1.137', 'B.1.70', 'B.1.1.148', 'B.1.190', 'B.1.1.289', 'B.1.249', 'B.1.160.28', 'B.1.1.161', 'B.1.1.170', 'B.1.1.429', 'B.1.22.1', 'B.1.389', 'B.1.527', 'B.1.1.153', 'B.1.258.11', 'B.1.258.21', 'B.1.1.44', 'B.1.160.8', 'B.1.177.87', 'C.40', 'B.1.397', 'B.1.413', 'B.1.1.187', 'B.1.1.229', 'B.1.1.406', 'B.1.177.33', 'C.18', 'B.1.177.17', 'B.1.1.307', 'B.1.1.311', 'B.1.177.7', 'B.1.36.39', 'B.1.231', 'B.1.1.240', 'B.1.177.16', 'B.1.177.57', 'B.1.549', 'B.1.177.84', 'C.38', 'B.1.177.24', 'B.1.468', 'B.1.177.54', 'B.1.177.12', 'B.1.177.69', 'B.1.160.7', 'B.1.177.6', 'B.1.1.46', 'B.1.408', 'B.1.177.19', 'B.1.1.130', 'B.1.258.5', 'B.1.1.241', 'B.1.177.9', 'Z.1', 'B.1.36.28', 'AP.1', 'B.1.362.2', 'B.1.258.7', 'B.1.221.1', 'B.1.177.56', 'B.1.177.20', 'B.1.177.48', 'B.1.1.216', 'B.5', 'B.1.1.283', 'B.1.1.226', 'B.6', 'B.1.282', 'B.1.1.48', 'B.1.617.1', 'C.36.3', 'B.1.36.16', 'B.1.1.51', 'B.1.351.3', 'AY.1', 'AY.74', 'AY.43', 'BA.1.18', 'B.1.178', 'B.1.273', 'B.1.6', 'B.1.225', 'B.1.1.409', 'B.1.1.142', 'B.1.78', 'B.1.1.362', 'B.1.94', 'B.1.8', 'B.1.1.58', 'B.1.1.464', 'B.1.1.47', 'B.1.9.4', 'B.38', 'B.4', 'B.1.9.5', 'B.1.350', 'B.1.96', 'B.1.1.185', 'B.42', 'B.1.1.56', 'B.1.1.301', 'B.1.1.305', 'B.1.218', 'B.1.1.221', 'B.1.1.269', 'B.1.1.33', 'B.1.416', 'B.1.597', 'N.2', 'B.1.219', 'AH.2', 'B.1.1.433', 'B.1.1.266', 'B.1.1.50', 'B.1.1.521', 'B.1.180', 'B.1.1.218', 'B.1.240', 'C.5', 'B.1.1.144', 'C.30', 'B.1.609', 'B.1.1.282', 'L.2', 'B.1.310', 'B.1.453', 'B.1.445', 'B.1.274', 'B.1.195', 'B.1.437', 'B.1.206', 'B.1.177.31', 'B.1.1.74', 'B.1.258.9', 'B.1.1.326', 'B.1.160.11', 'B.1.1.242', 'B.1.467', 'B.1.1.297', 'AN.1', 'B.1.398', 'B.1.1.243', 'B.1.329', 'B.1.36.9', 'B.1.177.26', 'B.1.36.2', 'B.1.177.8', 'B.1.1.57', 'B.1.1.111', 'B.1.1.459', 'B.1.381', 'B.1.1.34', 'B.1.429', 'B.41', 'B.43', 'B.1.371', 'B.1.595', 'B.30', 'B.1.108', 'A.3', 'B.1.38', 'B.1.1.61', 'B.1.450', 'B.1.314', 'B.1.211', 'B.1.37', 'B.1.110', 'B.1.320', 'B.4.4', 'B.1.369', 'B.1.319', 'B.1.162', 'B.1.446', 'B.1.166', 'B.1.370', 'B.1.264', 'B.1.473', 'B.1.378', 'B.1.518', 'B.1.586', 'B.1.333', 'B.1.302', 'B.1.507', 'B.6.6', 'B.1.169', 'B.1.494', 'B.1.323', 'B.1.360', 'B.1.119', 'B.53', 'B.1.1.306', 'B.4.1', 'B.1.110.1', 'B.1.210', 'B.1.564', 'B.1.36.8', 'B.1.145', 'B.1.14', 'B.20', 'B.19', 'A.4', 'B.1.384', 'B.1.112', 'B.1.426', 'B.1.552', 'B.1.301', 'B.1.199', 'A.6', 'B.46', 'B.4.5', 'B.1.1.416', 'B.1.23', 'B.1.338', 'B.1.1.16', 'B.1.503', 'B.1.247', 'B.1.324', 'B.1.276', 'B.1.103', 'B.1.139', 'B.1.110.3', 'B.1.293', 'B.1.143', 'B.1.187', 'B.1.528', 'B.1.268', 'B.1.390', 'N.1', 'B.1.479', 'B.1.415', 'B.1.435', 'B.1.442', 'B.1.1.231', 'B.1.523', 'B.1.441', 'B.1.260', 'B.1.1.319', 'B.1.1.487', 'B.1.179', 'B.1.1.175', 'B.1.515', 'B.1.485', 'B.1.311', 'B.1.164', 'D.2', 'B.6.3', 'B.1.170', 'C.17', 'B.1.336', 'B.1.1.26', 'B.1.124', 'B.1.163', 'B.1.1.113', 'B.1.241', 'B.1.558', 'B.1.334', 'B.1.113', 'B.1.434', 'B.4.6', 'B.1.495', 'B.1.520', 'B.1.9.2', 'B.1.509', 'B.1.483', 'B.1.578', 'B.51', 'B.1.521', 'B.13', 'A.19', 'B.1.281', 'B.1.1.59', 'B.1.1.28', 'AE.2', 'AE.4', 'B.1.1.136', 'B.1.471', 'B.1.313', 'B.1.1.133', 'D.3', 'B.1.330', 'B.1.3', 'B.1.576', 'B.1.304', 'B.1.294', 'B.1.401', 'B.1.1.172', 'B.1.1.132', 'B.1.508', 'B.1.159', 'B.1.493', 'B.1.428', 'B.1.106', 'B.1.425', 'B.1.1.437', 'B.1.315', 'B.1.1.181', 'B.1.498', 'B.1.501', 'B.1.366', 'B.1.1.192', 'B.37', 'B.1.189', 'B.1.341', 'B.1.400', 'B.1.267', 'B.1.245', 'B.1.399', 'B.1.405', 'B.1.239', 'B.1.556', 'B.1.1.380', 'B.1.126', 'B.1.1.128', 'B.1.1.93', 'B.1.1.290', 'B.1.328', 'B.1.594', 'B.1.1.158', 'B.1.554', 'B.1.1.116', 'B.1.306', 'B.1.1.244', 'B.1.1.180', 'B.1.595.1', 'B.1.502', 'B.1.340', 'B.1.182', 'B.1.318', 'B.1.1.225', 'B.1.1.228', 'B.1.1.376', 'B.1.448', 'B.1.1.135', 'B.1.565', 'B.1.305', 'B.1.349', 'B.1.385', 'A.12', 'B.1.1.335', 'B.1.359', 'B.1.1.67', 'A.11', 'B.1.1.404', 'B.1.1.441', 'B.1.280', 'B.1.451', 'B.1.1.298', 'B.1.1.63', 'B.1.1.263', 'B.1.1.526', 'B.1.595.3', 'B.1.574', 'B.1.582', 'B.1.265', 'B.1.423', 'B.1.436', 'B.1.605', 'B.1.428.1', 'B.1.606', 'B.1.1.237', 'B.1.181', 'B.1.335', 'B.1.1.98', 'B.1.443', 'B.1.547', 'B.1.184', 'B.1.36.10', 'B.1.1.337', 'B.1.325', 'B.1.580', 'B.1.382', 'B.1.595.4', 'B.1.116', 'B.1.487', 'B.1.1.186', 'A.25', 'B.1.1.431', 'B.1.36.38', 'B.1.188', 'B.1.452', 'B.1.496', 'B.1.1.230', 'B.1.298', 'B.1.1.265', 'B.1.403', 'B.1.287', 'B.1.1.328', 'B.1.612', 'B.1.570', 'B.1.1.77', 'C.14', 'B.1.205', 'C.13', 'C.33', 'C.32', 'C.4', 'C.25', 'B.1.1.381', 'B.1.1.110', 'B.1.396', 'B.1.1.251', 'C.31', 'B.1.337', 'C.23', 'B.1.1.432', 'B.1.377', 'B.1.234', 'B.1.1.512', 'B.1.596', 'B.1.513', 'B.1.1.205', 'B.1.1.423', 'B.1.355', 'B.1.492', 'B.1.602', 'B.1.1.207', 'B.1.444', 'B.1.284', 'B.1.566', 'B.1.589', 'B.1.590', 'B.1.1.447', 'B.1.546', 'B.1.544', 'B.1.577', 'B.1.1.320', 'B.1.567', 'B.1.361', 'B.1.421', 'B.1.500', 'B.1.1.440', 'B.1.1.356', 'B.1.1.210', 'B.1.1.258', 'B.1.1.367', 'B.1.1.118', 'B.1.1.402', 'B.1.1.268', 'B.1.572', 'B.1.1.463', 'B.1.1.139', 'B.1.309', 'B.1.1.239', 'B.1.1.514', 'B.1.599', 'B.1.432', 'B.1.561', 'B.1.1.329', 'B.1.36.27', 'B.1.36.18', 'B.1.499', 'B.1.194', 'B.1.571', 'B.1.424', 'B.1.1.513', 'B.1.504', 'B.1.516', 'B.1.362', 'B.1.1.453', 'B.1.1.434', 'B.1.600', 'B.1.541', 'B.1.1.392', 'AE.7', 'B.1.232', 'AE.6', 'B.1.1.370', 'B.1.533', 'B.1.428.2', 'B.1.213', 'B.1.402', 'L.4', 'B.1.1.344', 'B.1.1.273', 'B.1.1.169', 'B.1.596.1', 'B.1.1.203', 'N.7', 'B.1.1.403', 'B.1.505', 'B.1.1.121', 'B.1.36.26', 'B.1.1.176', 'B.1.1.424', 'B.1.1.27', 'B.1.591', 'B.1.497', 'N.6', 'N.4', 'C.26', 'B.1.134', 'B.1.409', 'B.1.568', 'B.1.476', 'B.1.517', 'B.1.375', 'B.1.1.342', 'B.1.511', 'B.1.1.177', 'B.1.233', 'B.1.354', 'B.1.1.303', 'B.1.36.24', 'B.1.1.398', 'B.1.581', 'B.1.1.368', 'B.1.1.467', 'B.1.404', 'B.1.1.316', 'B.1.535', 'B.1.551', 'B.1.1.322', 'B.1.433', 'B.1.478', 'B.1.543', 'B.1.1.486', 'B.1.1.352', 'B.1.258.23', 'B.1.1.312', 'B.1.550', 'B.1.601', 'Q.3', 'B.1.1.334', 'B.1.1.351', 'B.1.1.426', 'B.1.1.389', 'B.1.587', 'B.1.411', 'B.1.420', 'W.1', 'B.1.400.1', 'B.1.603', 'B.1.1.117', 'B.1.1.353', 'B.1.526', 'B.1.438.4', 'B.1.1.348', 'B.1.595.2', 'B.1.506', 'B.1.1.359', 'B.1.482', 'A.2.5', 'B.1.539', 'B.1.1.480', 'B.1.118', 'B.1.177.53', 'B.1.1.517', 'B.1.1.518', 'B.1.1.159', 'B.1.348', 'B.1.1.418', 'B.1.36.36', 'B.1.1.450', 'B.1.343', 'B.1.499.1', 'N.5', 'B.1.560', 'B.1.635', 'B.1.351.2', 'R.2', 'B.1.177.18', 'B.1.636', 'B.1.357', 'B.1.1.261', 'B.1.1.393', 'B.1.1.411', 'B.1.623', 'B.1.429.1', 'A.2.5.2', 'B.1.1.461', 'B.1.438.3', 'N.3', 'AE.8', 'B.1.564.1', 'B.1.585', 'B.1.563', 'B.1.438', 'B.1.1.294', 'B.1.177.46', 'B.1.192', 'B.1.1.200', 'B.1.569', 'B.1.530', 'B.1.1.54', 'B.1.1.318', 'C.37', 'B.1.637', 'B.1.279', 'B.1.363', 'B.1.545', 'B.1.36.22', 'B.1.438.1', 'B.1.36.7', 'P.1.10', 'B.1.243.1', 'B.1.1.345', 'B.1.277', 'B.36', 'B.1.115', 'C.6', 'B.1.1.354', 'B.1.631', 'B.1.604', 'B.1.1.315', 'B.1.1.413', 'B.1.1.340', 'B.1.625', 'XB', 'B.1.463', 'B.1.237', 'B.1.627', 'B.1.1.134', 'B.1.1.358', 'AT.1', 'N.8', 'B.1.1.152', 'B.1.524', 'B.1.177.75', 'B.1.1.374', 'B.1.1.99', 'B.1.517.1', 'A.7', 'B.1.369.1', 'B.1.84', 'P.1.17', 'P.1.12', 'P.1.13', 'P.1.1', 'B.1.137', 'A.24', 'Q.2', 'B.1.1.163', 'B.1.466.2', 'P.1.15', 'B.1.555', 'B.1.1.462', 'B.1.634', 'B.1.575.1', 'P.3', 'C.37.1', 'B.1.1.141', 'P.1.14', 'B.1.620', 'C.2.1', 'XH', 'A.2.5.1', 'AZ.3', 'A.29', 'B.1.1.525', 'XG', 'B.1.243.2', 'B.1.36.34', 'AY.65', 'AY.61', 'AY.122', 'AY.102', 'P.1.2', 'C.28', 'A.2.5.3', 'B.1.1.262', 'B.1.1.288', 'P.7', 'B.1.1.387', 'AY.16', 'B.1.1.272', 'B.1.621', 'C.36.3.1', 'B.1.1.523', 'AY.44', 'AY.28', 'AY.75', 'Q.8', 'AY.38', 'B.1.177.83', 'Q.6', 'AY.9', 'AY.20', 'B.1.1.515', 'AY.49', 'P.1.16', 'AY.110', 'AY.25', 'AY.37', 'AY.26', 'AY.4.8', 'AY.17', 'AY.100', 'AY.66', 'AY.55', 'AY.112', 'AY.52', 'B.1.626', 'AY.126', 'AY.9.2', 'AY.39', 'AY.116', 'AY.10', 'AY.25.1', 'AY.13', 'AY.62', 'AY.127', 'AY.3', 'AY.119', 'N.9', 'AY.48', 'AY.99', 'AY.70', 'B.1.1.204', 'B.1.630', 'AY.71', 'AY.47', 'B.1.1.122', 'AY.67', 'AY.129', 'B.1.1.407', 'AY.72', 'AY.56', 'B.1.638', 'C.36.1', 'AY.4', 'AY.54', 'AY.120.2', 'B.1.1.373', 'AY.92', 'AY.36.1', 'A.23', 'B.1.393', 'B.1.214', 'AY.120', 'P.1.17.1', 'B.1.36.17', 'P.1.7.1', 'AY.60', 'AY.19', 'B.1.289', 'AY.2', 'AY.103', 'XR', 'AY.98.1', 'AY.39.1', 'AY.64', 'AY.117', 'P.1.7', 'B.1.177.80', 'AY.114', 'AY.77', 'AY.46', 'AY.119.2', 'AY.107', 'P.6', 'B.1.617.3', 'B.1.621.1', 'B.1.1.115', 'AY.35', 'A.2.4', 'B.1.161', 'AY.24', 'AY.120.1', 'AY.105', 'AY.5', 'AY.3.1', 'AY.121', 'AY.46.4', 'AY.59', 'AY.50', 'AY.34.1.1', 'AY.106', 'AY.86', 'AY.16.1', 'AY.40', 'AY.118', 'B.1.593', 'AY.15', 'AY.81', 'AY.76', 'AY.51', 'AY.53', 'AY.58', 'B.1.415.1', 'AY.124', 'AY.5.4', 'B.1.1.371', 'B.1.1.397', 'AY.33', 'AY.25.3', 'AY.45', 'AY.87', 'AY.46.2', 'AY.125', 'AY.46.6', 'B.1.632', 'AY.46.1', 'AY.116.1', 'AY.22', 'AY.78', 'AY.91', 'AY.98', 'AY.5.3', 'AY.109', 'AY.113', 'AY.95', 'AY.124.1.1', 'AY.42', 'B.1.575.2', 'AY.121.1', 'AY.32', 'AY.7.2', 'AY.80', 'AY.34', 'AY.108', 'AY.4.7', 'AY.68', 'AY.112.1', 'B.1.1.382', 'AY.103.2', 'BB.2', 'AY.83', 'AY.88', 'AY.3.4', 'AY.73', 'AY.4.2.3', 'AY.39.2', 'AY.36', 'AY.90', 'AY.128', 'AY.69', 'B.1.538', 'AY.27', 'B.1.621.2', 'AY.25.2', 'AY.42.1', 'AY.6', 'AY.93', 'AY.82', 'B.1.1.82', 'AY.111', 'AY.5.5', 'AY.57', 'AY.3.3', 'AY.4.9', 'AY.123', 'AY.99.2', 'AY.103.1', 'AY.46.3', 'AY.34.2', 'AY.122.2', 'AY.122.1', 'AY.4.6', 'AY.102.1', 'AY.7', 'AY.132', 'AY.23.1', 'AY.24.1', 'AY.94', 'AY.25.1.1', 'AY.5.1', 'AY.79', 'AY.101', 'B.1.629', 'AY.11', 'AY.5.6', 'AY.122.4', 'AY.39.3', 'AY.30', 'B.1.579', 'B.1.422', 'B.1.407', 'B.1.395', 'B.1.1.72', 'B.1.607', 'B.1.1.422', 'AY.41', 'B.1.1.331', 'B.1.292', 'AY.43.8', 'AY.7.1', 'AY.5.2', 'P.1.12.1', 'AY.39.1.4', 'B.1.639', 'B.1.431', 'B.1.69', 'B.1.9.1', 'B.1.168', 'B.1.40', 'B.1.1.160', 'B.1.1.14', 'B.1.418', 'B.1.252', 'B.1.1.147', 'B.1.1.349', 'W.3', 'B.1.36.12', 'B.1.1.138', 'B.1.1.309', 'D.4', 'B.1.235', 'B.1.36.23', 'B.1.1.286', 'B.1.1.253', 'B.1.258.2', 'B.1.610', 'B.1.1.280', 'B.1.1.255', 'B.1.258.3', 'B.1.1.227', 'B.4.8', 'B.1.532', 'B.1.1.365', 'C.2', 'B.1.1.107', 'B.1.1.410', 'B.1.1.234', 'B.1.1.347', 'B.1.238', 'B.1.1.408', 'B.1.1.421', 'B.1.77', 'B.1.1.12', 'B.1.1.97', 'B.1.1.299', 'B.1.177.65', 'B.1.1.120', 'B.1.1.15', 'C.3', 'B.1.1.145', 'B.1.1.125', 'B.1.1.300', 'B.1.1.363', 'B.1.1.256', 'AK.1', 'B.1.1.302', 'B.1.173', 'A.26', 'B.1.110.2', 'B.1.1.162', 'B.1.1.209', 'B.1.1.310', 'B.1.177.10', 'B.1.1.196', 'B.1.1.165', 'B.1.177.3', 'B.1.1.364', 'B.1.1.425', 'AE.1', 'B.1.1.361', 'B.1.1.109', 'B.1.1.92', 'B.1.1.155', 'B.1.1.123', 'M.3', 'B.1.469', 'S.1', 'B.1.251', 'B.1.1.95', 'B.1.1.55', 'B.1.372', 'B.1.1.149', 'B.1.1.86', 'B.1.350.1', 'B.1.1.38', 'AE.5', 'B.1.1.201', 'B.1.1.325', 'B.1.510', 'C.1', 'B.18', 'B.1.1.154', 'B.1.1.308', 'B.1.177.2', 'B.1.177.55', 'B.1.258.4', 'G.1', 'B.1.160.23', 'B.1.177.11', 'B.1.1.236', 'B.1.177.5', 'B.1.258.6', 'B.1.177.38', 'B.1.460', 'B.1.258.10', 'B.1.240.2', 'B.1.1.341', 'B.1.177.58', 'B.1.258.12', 'B.1.362.1', 'AS.2', 'B.1.177.14', 'B.1.177.50', 'B.1.1.193', 'B.1.1.500', 'B.1.177.35', 'B.1.428.3', 'B.1.177.21', 'B.1.1.375', 'AC.1', 'AA.4', 'AA.6', 'B.1.160.33', 'AS.1', 'B.1.177.30', 'B.1.177.63', 'AA.1', 'B.1.177.41', 'AH.3', 'V.1', 'B.1.258.22', 'B.1.36.21', 'AA.3', 'AA.7', 'B.1.221.2', 'B.1.177.74', 'B.1.36.33', 'AA.2', 'B.1.1.355', 'B.1.562', 'B.1.1.182', 'B.1.177.47', 'B.1.462', 'B.1.177.73', 'B.1.459', 'B.1.177.67', 'B.1.1.417', 'B.1.531', 'U.1', 'B.1.177.45', 'D.5', 'AM.1', 'AA.5', 'B.1.177.59', 'AA.8', 'B.1.208', 'B.1.258.20', 'B.1.36.20', 'B.1.36.37', 'XA', 'AF.1', 'C.36.2', 'B.1.160.20', 'B.1.177.40', 'B.1.160.19', 'B.1.177.89', 'B.1.1.428', 'B.1.1.327', 'B.1.1.430', 'B.1.1.202', 'C.20', 'AH.1', 'B.1.1.385', 'B.1.1.458', 'B.1.1.336', 'AB.1', 'B.1.264.1', 'AY.4.5', 'AY.99.1', 'P.1.9', 'AY.8', 'B.1.351.5', 'AY.84', 'AY.120.2.1', 'B.1.380', 'AY.25.1.2', 'AY.39.1.1', 'P.1.4', 'AY.75.2', 'AY.39.1.2', 'AY.127.1', 'AZ.2', 'AY.4.4', 'AY.3.2', 'B.1.1.524', 'AY.134', 'B.1.618', 'N.10', 'AY.4.2', 'AY.112.3', 'AY.43.4', 'A.9', 'AY.4.1', 'C.12', 'B.1.1.412', 'B.1.1.506', 'C.1.2', 'B.6.8', 'AY.104', 'B.1.1.448', 'AZ.5', 'AY.131', 'B.50', 'AY.122.5', 'B.1.344', 'AY.43.1', 'AY.46.5', 'AY.5.7', 'AY.33.1', 'AY.43.3', 'AY.4.13', 'AY.33.2', 'B.1.617', 'AY.20.1', 'AY.122.6', 'AY.112.2', 'AY.4.11', 'B.1.640.1', 'AY.4.2.1', 'AY.4.3', 'AY.4.2.4', 'AY.119.1', 'B.1.637.1', 'B.1.212', 'AY.31', 'AY.4.2.2', 'AY.39.1.3', 'AY.26.1', 'AY.127.3', 'AU.3', 'AZ.1', 'B.1.537', 'AY.4.10', 'B.1.1.528', 'BA.1.10', 'BA.1.13', 'AY.4.17', 'AY.124.1', 'B.1.640.2', 'BA.1.19', 'BA.1.21', 'B.1.351.1', 'AW.1', 'BA.1.17', 'AY.18', 'BA.1.1.14', 'BA.1.1.8', 'BA.1.1.10', 'AY.4.15', 'B.1.263', 'B.49', 'BA.1.1.18', 'AY.43.5', 'B.1.36.19', 'BA.1.15.1', 'BA.1.1.6', 'BA.1.15.2', 'BA.1.1.17', 'BA.1.9', 'BA.1.1.16', 'BA.1.8', 'BA.1.1.13', 'BA.1.1.15', 'BA.2.9', 'B.1.548', 'BA.1.14', 'BA.1.16', 'AY.46.6.1', 'AY.133', 'BA.1.6', 'AD.1', 'BA.1.14.2', 'BA.1.17.1', 'BA.1.1.4', 'BA.1.14.1', 'BA.1.1.11', 'B.1.540', 'BA.1.1.12', 'BA.1.7', 'BA.1.1.7', 'XS', 'BA.1.21.1', 'BA.1.3', 'BA.1.12', 'BA.2.12', 'BA.3', 'AY.43.9', 'AY.43.7', 'BA.2.25', 'BA.2.23', 'B.1.1.157', 'BA.2.3.3', 'BA.2.27', 'BA.2.7', 'BA.2.4', 'B.1.1.40', 'B.1.1.383', 'BA.1.5', 'B.1.1.62', 'B.1.1.507', 'B.1.1.53', 'AU.1', 'BA.2.9.2', 'BA.2.5', 'BA.1.13.1', 'BA.2.9.1', 'BA.2.32', 'BA.2.1', 'BA.2.6', 'BA.2.2', 'C.21', 'B.1.1.388', 'B.1.1.332', 'P.1.8', 'BA.2.20', 'BA.2.21', 'BA.2.12.1', 'BA.2.18', 'BA.2.26', 'BA.3.1', 'BA.2.8', 'BA.2.3.4', 'BA.2.31', 'BA.2.14', 'BA.2.17', 'XE', 'AY.9.2.2', 'BA.2.15', 'BA.2.34', 'BA.1.1.3', 'C.8', 'BA.2.30', 'BA.2.19', 'BA.2.13', 'BA.2.22', 'BA.2.11', 'BA.2.28', 'BA.4', 'BA.5', 'BA.2.16', 'M.2', 'B.1.177.64', 'B.1.326', 'AZ.2.1', 'B.1.1.366', 'U.3', 'B.1.177.49', 'B.1.1.396', 'B.1.221.3', 'AG.1', 'B.1.160.13', 'AK.2', 'B.1.1.224', 'B.1.1.442', 'W.2', 'B.1.1.338', 'B.1.1.129', 'B.1.221.4', 'B.1.258.24', 'B.1.160.25', 'Q.7', 'B.47', 'B.1.387', 'AQ.1', 'B.1.177.42', 'B.1.1.333', 'AV.1', 'AZ.4', 'AY.21', 'B.1.1.446', 'B.1.633', 'AY.91.1', 'B.1.1.449', 'B.1.177.39', 'B.1.489', 'AM.3', 'AY.102.2', 'B.1.619.1', 'AY.63', 'AY.4.14', 'AY.23.2', 'AY.4.16', 'AY.43.2', 'AY.4.12', 'B.1.1.88', 'AY.43.6', 'AY.4.2.5', 'AY.123.1', 'BA.1.4', 'AY.125.1', 'BA.1.16.2', 'BA.1.16.1', 'AY.122.3', 'BA.1.22', 'B.1.177.70', 'BA.2.33', 'BA.2.25.1']
    #print(avail_lineages)
    #accessions = dc.getAccessionsByLineage(avail_lineages[0])[:7] #### take 7 accessions of B.1.1 lineage 
    accessions = ['OA977655', 'OD917797', 'OA971744', 'OU081140', 'FR995488', 'MT451035', 'MW705321']
    #print(accessions)
    
    dict = {}        
    directory  = './data/raw'  # <--- substitute by your local path
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f): 
            sequence = extractFasta(f)
            kmers_list = kmers(4, 1, sequence)
            print(kmers_list)
            dict[filename] = kmers_list
    #print(dict)
    print(database1_build(dict))
    reps_of_lineage = lineage_reps_from_its_accessions(dict,13)
            
    database1={}    #global dictionary for database 1
    database1.update({avail_lineages[0]:reps_of_lineage})
    print(database1)
