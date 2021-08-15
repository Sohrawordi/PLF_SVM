# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 21:35:46 2021

@author: Md. Sohrawordi
"""
from flask import Flask, render_template,request,send_file


def amino_acid_composition(dataset):
    import pandas as pd
    #A list of 21 Amino Acid
    amino_acid="ACDEFGHIKLMNPQRSTWYVX"
    colls=[]  
    #sample = pd.read_csv(dataset)
    
    for i in range (len(amino_acid)):
        colls.append("Amino_acid_"+amino_acid[i])
    #print(colls)
    df=pd.DataFrame(columns=colls)
    for i in range(len(dataset['Sequence'])):
        #print(i)
        row={}
        for acid in amino_acid:
            #print(seq.count(acid))
            row["Amino_acid_"+acid]=dataset['Sequence'][i].count(acid)/len(dataset['Sequence'])
        df=df.append(row,ignore_index=True)

    



    return df


def binary_data(positive_sample):
    import pandas as pd
    amino_acid="ACDEFGHIKLMNPQRSTWYVX"
    binary_code={}
    colls=[]    
    for i in range (len(amino_acid)*len(positive_sample['Sequence'][0])):
        colls.append("Binary"+str(i))

    for i in amino_acid:
        a=[]
        for j in amino_acid:
            if i==j:
                a.append(1)
            else:
                a.append(0)  
        binary_code[i]=a 
        
        
        
    pdf=pd.DataFrame(columns=colls)
    for k in range(len(positive_sample['Sequence'])):
        c=0 
        row_value={}
        for i in positive_sample['Sequence'][k]:
            for j in binary_code[i]:
                row_value["Binary"+str(c)]=j
                c=c+1
        ##row_value["Class"]=1        
        pdf=pdf.append(row_value,ignore_index=True) 
    
    return pdf






def count(seq,a,b):
    fp=seq.find(a)
    fpb=seq.find(b)
    
    
    s=0
    if fp>=0 and fpb>=0:
        for j in range(len(seq)):
            if b==seq[j]:
                s=s+abs(j-fp)
            
        #ar.append(s)
        return s
    else:
        return -1
        #ar.append(-1)
    #return ar




def prim_rprim(dataset):
    import pandas as pd
    
    amino_acid="ACDEFGHIKLMNPQRSTWYVX"
    
    
    colls=[]       
    for x in amino_acid:
        for y in amino_acid:
            colls.append("prim_"+x+"-->"+y)
   
    pdf=pd.DataFrame(columns=colls)    
    
    
    for x in amino_acid:
        for y in amino_acid:
            colls.append("rv_prim_"+x+"-->"+y)            
            
    colls=[]               
    rdf=pd.DataFrame(columns=colls)
    
    
    #calculation of max length between resudues
    seq_length=len(dataset['Sequence'][0])
    max_distance=0
    for i in range(seq_length):
        max_distance=max_distance+i

    #print(max_distance)
    
    
    
    
    
    
    
    
    
    
    
    for i in range(len(dataset.index)):
        #print(dataset['Sequence'][i])
        seq=dataset['Sequence'][i]
        #seq="ADCDD"
        
        
        
        pr_row={}
        
        
        # for position relative incidence matrix
        for x in amino_acid:
            for y in amino_acid:
                li=count(seq,x,y)
                #print(x,"....>",y,li)
                pr_row["prim_"+x+"-->"+y]=li/max_distance
                
            #print()
            #print()
          
            
        # for reverse position relative incidence matrix    
        rv_seq=seq[::-1] 
        r_row={}
        for x in amino_acid:
            for y in amino_acid:
                li=count(rv_seq,x,y)
                #print(x,"....>",y,li)
                r_row["rv_prim_"+x+"-->"+y]=li /max_distance   
            
                
        #print(dataset['Sequence'][i])
 
        pdf=pdf.append(pr_row,ignore_index=True)
        rdf=rdf.append(r_row,ignore_index=True)        
        
    return pdf,rdf




# split each sequence into characters and create a dataframe with a column for each character

def propensity_matrix_generation(dataset,window_size):

    amino_acid="ACDEFGHIKLMNPQRSTWYVX"
    strim_len=int(window_size/2)
    
    p_matrix={}
    #calculate total element of a position in dataset
    n=len(dataset.index)

    for i in range (window_size):
        ar_list=[]
        if i!=strim_len:
            for j in range(len(dataset.index)):
                #print(dataset['Sequence'][j][i])
                ar_list.append(dataset['Sequence'][j][i])
                
            for acid in amino_acid:
                #print(acid,ar_list.count(acid))
                p_matrix[acid+"_"+str(i-strim_len)]=ar_list.count(acid)/n
                
        
        

            
            
        
        
    
    
    
    return p_matrix
    


def propensity_feature(dataset,window_size):
    import pickle
    matrix = pickle.load(open('propensity_matrix.pkl','rb'))

    import pandas as pd
    #amino_acid="ACDEFGHIKLMNPQRSTWYVX"
    strim_len=int(window_size/2)
    colls=[]
    #generate dataframe columns without central resudue
    for i in range (window_size):
        if i !=strim_len:
            colls.append("propensity_"+str(i-strim_len))
    #print(colls)
    df=pd.DataFrame(columns=colls)
    
    
    
    for i in range(len(dataset['Sequence'])):
        #print(dataset['Sequence'][i])
        seq=dataset['Sequence'][i]
        row={}
        
        #create a dataframe row without central resudue for each sample sequence
        for j in range(len(seq)):
            if j !=strim_len:
                row["propensity_"+str(j-strim_len)]=matrix[seq[j]+"_"+str(j-strim_len)]
       
        #add a row in dataframe
        df=df.append(row,ignore_index=True)


    return df


def featureset(data):
    import pickle
    import pandas as pd
    aac=amino_acid_composition(data)
    be=binary_data(data)
    prim,rprim=prim_rprim(data)
    psaap=propensity_feature(data,21)
    
    


    hybrid_sample=aac.copy()
    for col in be.columns: 
        hybrid_sample[col]=be[col] 
    for col in psaap.columns:
        hybrid_sample[col]=psaap[col] 
    for col in prim.columns: 
        hybrid_sample[col]=prim[col]  

    for col in rprim.columns: 
        hybrid_sample[col]=rprim[col]
        
    
    
    
    cols_index=pickle.load(open('optimal_feature_set_index.pkl','rb'))

    optimal_data=pd.DataFrame()
    for index in cols_index:
        optimal_data[index]=hybrid_sample[index]       
      
        
    #print(optimal_data)  
    return optimal_data





def prediction(seq,window):
    up=int(window/2)

    listk=[i for i, letter in enumerate(seq) if letter == 'K']
    import pandas as pd
    df=pd.DataFrame(columns=["Sequence","Position"])
    #print(listk)
    for i in listk:
        if i<up:
            #print(seq[0:i+up+1],i)
            raw={}
            tem=""
            for j in range(window-len(seq[0:i+up+1])):
                tem=tem+'X'  
            newseq=tem+seq[0:i+up+1]
            raw["Sequence"]=newseq
            raw["Position"]=i+1
            df=df.append(raw,ignore_index=True)
            #print(newseq,len(newseq),i)

        
        else:
            raw={}
            newseq=seq[i-up:i+up+1]
            #print(newseq,len(newseq),i)
            raw["Sequence"]=newseq
            raw["Position"]=i+1
            df=df.append(raw,ignore_index=True)
            
    x=featureset(df)
    train=pd.read_csv('all_samples_optimal_features.csv') 
    label=train.pop("Class")
    from sklearn.svm import SVC
    classifier =SVC(kernel='rbf', C=4,gamma=1,probability=True) 
    model=classifier.fit(train, label)
    pred=model.predict(x)
    df['Class']=pred
        
    
    return df



app = Flask(__name__)

@app.route("/")
def home():

    return render_template('index.html')
@app.route('/about')
def about_me():
    return render_template('about.html')

@app.route('/download')
def download_file():
    file='raw_data.pdf'
    return send_file(file,as_attachment=True)



@app.route('/help')
def help():
    return render_template('help.html')


@app.route('/predict',methods=['POST']) 
def predict():
    if request.method=="POST":
        
        
        seq=request.form["protein"]
        listk=[i for i, letter in enumerate(seq) if letter == 'K']
        if len(seq)>20 and len(listk)>0:
            #if seq[len(seq)-1]=='K':
            #seq=seq+"A"
            seq=seq+"XXXXXXXXXXXXXXXXXXXXXXXXX"
            result=prediction(seq,21)
            result['Class'].replace({0: "Non-Formylated Lysine", 1: "Formylated Lysine"}, inplace=True)

            return render_template('show.html',n=result) 

        else:
            return "The length of the given protein may be less than 21 <br> You have entered blank space<br>   The protein does not contain any Lysine residue<br> <b>Please try again with another larger protein</b>"

        
        
        






if __name__ == "__main__": 
    app.run(debug=True) 