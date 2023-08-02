# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

Notes
- 230731 Some of the parsed files include the records whose Active Substances is np.nan
    - 2014q4_0, 2017q2_0, 2019q1_2, 2018q4_2

@author: tadahaya
"""

import pandas as pd
import xml.etree.ElementTree as ET
import os
from pathlib import Path
from tqdm import tqdm

FORMAT = ["2014Q3-","2012Q4-2014Q2"]

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

# abstract class
class XMLoader():
    def __init__(self):
        self.__load = Loader2014Q3_()

    def get_format(self):
        return FORMAT

    def load_xml(self, url, to_pickle=False, to_csv=False, fileout="", sep="\t"):
        """
        a module for handling the large xml files of FAERS data
        take care the date of data format, which affects loading algorithm

        1) download zip files from FDA (https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)
        2) thaw a file into xml
        3) add the path to the below form
        4) run
        
        Parameters
        ----------
        url: str
            a path for xml file from FAERS
            
        to_pickle: boolean
            whether result is exported as pickle
            
        to_csv: boolean
            whether result is exported as csv
                    
        """
        return self.__load.load_xml(
            url, to_pickle=to_pickle, to_csv=to_csv, sep=sep, fileout=fileout
            )


    def to_2014Q3_(self):
        """ load xml data of 2014Q3- """
        self.__load = Loader2014Q3_()


    def to_2012Q4_2014Q2(self):
        """ load xml data of 2012Q4-2014Q2 """
        self.__load = Loader2012Q4_2014Q2()


# concrete class
class Loader2014Q3_():
    def __init__(self):
        pass

    def load_xml(self, url, to_pickle=False, to_csv=False, fileout="", sep="\t"):
        """
        a module for handling the large xml files of FAERS data
        it can be applicable to 2013-2018 data files so far
        
        1) download zip files from FDA (https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)
        2) thaw a file into xml
        3) add the path to the below form
        4) run
        
        Parameters
        ----------
        url: str
            a path for xml file from FAERS
            
        to_pickle: boolean
            whether result is exported as pickle
            
        to_csv: boolean
            whether result is exported as csv

        """
        if len(fileout)==0:            
            if to_pickle:
                fileout = os.path.dirname(url) + SEP + os.path.basename(url).replace('xml','pkl')
            elif to_csv:
                fileout = os.path.dirname(url) + SEP + os.path.basename(url).replace('xml','txt')   

        #file loading
        tree = ET.parse(url) #consume ~2GB
        root = tree.getroot()
        #root
        ##tag=inchicsr, attrib={lang:en}
        
        #2nd; saftyreport for each case
        #tag_2nd = [v.tag for v in root]
        #att_2nd = [v.attrib for v in root]
        saftyreports = [v for v in root] #the 1st one is the header
        del saftyreports[0]
        
        #case ID
        caseID = [v.find('safetyreportid').text for v in saftyreports]
        
        #event country
        event_country = []
        ap = event_country.append
        for v in saftyreports:
            ele = v.find('occurcountry')
            if ele is None:
                ap('')
            else:
                ap(ele.text)
        del ele,v
        
        #event date
        event_date = []
        ap = event_date.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('summary')
            if c is None:
                ap('')
            else:
                ap(c[0].text.replace('CASE EVENT DATE: ',''))
        del v,b,c
        
        #age
        age = []
        ap = age.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('patientonsetage')
            if c is None:
                ap('')
            else:
                ap(c.text)
        del v,b,c
        
        #sex
        sex = []
        ap = sex.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('patientsex')
            if c is None:
                ap('')
            else:
                if str(c.text)=='1':
                    ap('Male')
                else:
                    ap('Female')
        del v,b,c
        
        #adverse reactions
        reactions = []
        ap = reactions.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.findall('reaction')
            reaction = ""
            for w in c:
                d = w.find('reactionmeddrapt')
                if d is not None: # then, d must not be empty
                    if len(reaction)==0:
                        reaction += d.text
                    else:
                        reaction += ";" + d.text
            ap(reaction)
        del v,b,c
        reactions = [v.lower() for v in reactions]

        #suspected drugs
        drugs = []
        ap = drugs.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.findall('drug')
            drug = [] # drugs in a patient
            for w in c:
                d = w.find('activesubstance')
                if d is not None:
                    drug += d[0].text.split("\\")
            drug = list(set(drug))
            if len(drug)==0:
                ap("")
            elif len(drug)==1:
                ap(drug[0])
            else:
                drug2 = drug[0]
                for w in drug[1:]:
                    drug2 += ";" + w
                ap(drug2)
        del v, b, c
        drugs = [v.lower() for v in drugs]
        
        # qualification
        qualifications = []
        ap = qualifications.append
        candi_list = [
            'Physician','Pharmacist','Other Health Professional','Lawyer','Consumer or non-health proffesional'
            ]
        for v in saftyreports:
            b = v.find('primarysource')
            c = b.find('qualification')
            if c is None:
                ap('')
            else:
                ap(candi_list[int(c.text)-1])
        del v,b,c

        # export data
        data = pd.DataFrame(
            {'Case ID':caseID,
             'Active Substances':drugs,
             'Reactions':reactions,
             'Sex':sex,
             'Event Date':event_date,
             'Event Country':event_country,
             'Patient Age':age,
             'Qualification':qualifications,}
             )
        del tree,root,event_country,event_date,age,sex,caseID,drugs,saftyreports,reactions,qualifications
        data = data.reset_index(drop=True) # modify 230731, to keep case ID
        if to_pickle:
            data.to_pickle(fileout)
        elif to_csv:
            data.to_csv(fileout, index=False, sep=sep)
        return data


# concrete class
class Loader2012Q4_2014Q2():
    def __init__(self):
        pass

    def load_xml(self, url, to_pickle=False, to_csv=False, fileout="", sep="\t"):
        """
        code for handling the large xml files of FAERS data
        This can be applicable to 2012Q4-2014Q2 data files so far

        """
        if len(fileout)==0:            
            if to_pickle:
                fileout = os.path.dirname(url) + SEP + os.path.basename(url).replace('xml','pkl')
            elif to_csv:
                fileout = os.path.dirname(url) + SEP + os.path.basename(url).replace('xml','txt')    
        
        #file loading
        tree = ET.parse(url) #consume ~2GB
        root = tree.getroot()
        #root
        ##tag=inchicsr, attrib={lang:en}
        
        #2nd; saftyreport for each case
        #tag_2nd = [v.tag for v in root]
        #att_2nd = [v.attrib for v in root]
        saftyreports = [v for v in root] #the 1st one is the header
        del saftyreports[0]
        
        #case ID
        caseID = [v.find('safetyreportid').text for v in saftyreports]
        
        #event country
        event_country = []
        ap = event_country.append
        for v in saftyreports:
            ele = v.find('occurcountry')
            if ele is None:
                ap('')
            else:
                ap(ele.text)
        del ele,v
        
        # event date
        # 230801, summary is empty. maybe the query is wrong
        event_date = []
        ap = event_date.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('summary')
            if c is None:
                ap('')
            else:
                ap(c[0].text.replace('CASE EVENT DATE: ',''))
        del v,b,c
        
        #age
        age = []
        ap = age.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('patientonsetage')
            if c is None:
                ap('')
            else:
                ap(c.text)
        del v,b,c
        
        #sex
        sex = []
        ap = sex.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.find('patientsex')
            if c is None:
                ap('')
            else:
                if str(c.text)=='1':
                    ap('Male')
                else:
                    ap('Female')
        del v,b,c
        
        #adverse reactions
        reactions = []
        ap = reactions.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.findall('reaction')
            reaction = ""
            for w in c:
                d = w.find('reactionmeddrapt')
                if d is not None:
                    if len(reaction)==0:
                        reaction += d.text
                    else:
                        reaction += ";" + d.text
            ap(reaction)
        del v,b,c
        reactions = [v.lower() for v in reactions]


        # suspected drugs
        # 230801, drug is empty. maybe the query (activesubstance) is wrong
        drugs = []
        ap = drugs.append
        for v in saftyreports:
            b = v.find('patient')
            c = b.findall('drug')
            drug = [] # drugs in a patient
            for w in c:
                d = w.find('activesubstance')
                if d is not None:
                    drug += d[0].text.split("\\")
            drug = list(set(drug))
            if len(drug)==0:
                ap("")
            elif len(drug)==1:
                ap(drug[0])
            else:
                drug2 = drug[0]
                for w in drug[1:]:
                    drug2 += ";" + w
                ap(drug2)
        del v,b,c
        drugs = [v.lower() for v in drugs]
        

        #export data
        data = pd.DataFrame(
            {'Case ID':caseID,
             'Active Substances':drugs,
             'Reactions':reactions,
             'Sex':sex,
             'Event Date':event_date,
             'Event Country':event_country,
             'Patient Age':age}
             )
        del tree,root,event_country,event_date,age,sex,caseID,drugs,saftyreports,reactions
        data = data.reset_index(drop=True) # modify 230731, to keep case ID
        if to_pickle:
            data.to_pickle(fileout)
        elif to_csv:
            data.to_csv(fileout, index=False, sep=sep)
        return data