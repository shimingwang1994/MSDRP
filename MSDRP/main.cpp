#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

vector < vector <int> > association;
vector < vector <float> > similarity_drug;
vector < vector <float> > similarity_cellline;
string drug_name[251];
string drug_name_187[187];
vector <string> drug_name_target;
float similarity_drug_187[187][187];
vector < vector <float> > similarity_drug_target;
string cellline_name[990];
string cellline_name_962[962];
vector < vector <float> > similarity_cellline_962;

typedef struct
{
    vector <int> points;
    int id;
}cluster;
vector <cluster> clusters_drug;//the clusters of drug
vector <cluster> clusters_cellline;//the clusters of cellline
int cluster_update_drug[251];
int cluster_update_cellline[990];

typedef struct
{
    float sim;
    int id;
}St;

void input()
{
    ifstream in;
    in.open("input data\\drug-cellline_association.txt");
    for(int i=0;i<251;i++)
    {
        vector <int> temp;
        for(int j=0;j<990;j++)
        {
            int a;
            in>>a;
            temp.push_back(a);
        }
        association.push_back(temp);
        temp.clear();
    }
    in.close();

    ifstream in1;
    in1.open("input data\\drug_name.txt");
    int num1=0;
    while(in1.good())
    {
        getline(in1,drug_name[num1]);
        num1++;
    }
    in1.close();

    ifstream in2;
    in2.open("input data\\drug_name-187.txt");
    int num2=0;
    while(in2.good())
    {
        getline(in2,drug_name_187[num2]);
        num2++;
    }
    in2.close();

    ifstream in22;
    in22.open("input data\\drug_name-with-target.txt");
    string nnn;
    while(in22.good())
    {
        getline(in22,nnn);
        drug_name_target.push_back(nnn);
    }
    in22.close();

    ifstream in3;
    in3.open("input data\\similarity_drug-187.txt");
    for(int i=0;i<187;i++)
        for(int j=0;j<187;j++)
            in3>>similarity_drug_187[i][j];
    in3.close();

    ifstream in33;
    in33.open("input data\\drug-similarity-target-network.txt");
    for(int i=0;i<250;i++)
    {
        vector <float> temppppp;
        for(int j=0;j<250;j++)
        {
            float a;
            in33>>a;
            temppppp.push_back(a);
        }
        similarity_drug_target.push_back(temppppp);
        temppppp.clear();
    }
    in33.close();

    ifstream in4;
    in4.open("input data\\cell line-ID.txt");
    int num3=0;
    while(in4.good())
    {
        getline(in4,cellline_name[num3]);
        num3++;
    }
    in4.close();

    ifstream in5;
    in5.open("input data\\cell line-ID-962.txt");
    int num4=0;
    while(in5.good())
    {
        getline(in5,cellline_name_962[num4]);
        num4++;
    }
    in5.close();

    ifstream in6;
    in6.open("input data\\similarity_cellline-962.txt");
    for(int i=0;i<962;i++)
    {
        vector <float> temp1;
        for(int j=0;j<962;j++)
        {
            float t;
            in6>>t;
            temp1.push_back(t);
        }
        similarity_cellline_962.push_back(temp1);
        temp1.clear();
    }
    in6.close();
}

void cal_similarity_drug()
{
    float a=0;
    for(int i=0;i<251;i++)
        for(int j=0;j<990;j++)
            a=a+association[i][j]*association[i][j];
    float par=a/251;
    for(int i=0;i<251;i++)
    {
        vector <float> ttt;
        for(int j=0;j<251;j++)
        {
            float sumand=0;
            for(int k=0;k<990;k++)
            {
            if((association[i][k]==1)&&(association[j][k]==1))
                sumand++;
            }
            float par1=1/(sumand+1);
            float sum=0;
            for(int k=0;k<990;k++)
            {
                sum=sum+(association[i][k]-association[j][k])*(association[i][k]-association[j][k]);
            }
            ttt.push_back(exp(-(par1/par)*sum));
        }
        similarity_drug.push_back(ttt);
        ttt.clear();
    }
    for(int i=0;i<251;i++)
        for(int j=0;j<251;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<187;c++)
            if(drug_name_187[c]==drug_name[i])
                x=c;
        for(int d=0;d<187;d++)
            if(drug_name_187[d]==drug_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_drug_187[x][y]>similarity_drug[i][j]))
        {
            similarity_drug[i][j]=similarity_drug_187[x][y];
        }
    }

    for(int i=0;i<251;i++)
        for(int j=0;j<251;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<250;c++)
            if(drug_name_target[c]==drug_name[i])
                x=c;
        for(int d=0;d<250;d++)
            if(drug_name_target[d]==drug_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_drug_target[x][y]>similarity_drug[i][j]))
        {
            similarity_drug[i][j]=similarity_drug_target[x][y];
        }
    }
}

void cal_similarity_cellline()
{
    float a=0;
    for(int i=0;i<251;i++)
        for(int j=0;j<990;j++)
            a=a+association[i][j]*association[i][j];
    float par=a/990;
    for(int i=0;i<990;i++)
    {
        vector <float> temp2;
        for(int j=0;j<990;j++)
        {
            float sumand=0;
            for(int k=0;k<251;k++)
            {
            if((association[k][i]==1)&&(association[k][j]==1))
                sumand++;
            }
            float par1=1/(sumand+1);
            float sum=0;
            for(int k=0;k<251;k++)
            {
                sum=sum+(association[k][i]-association[k][j])*(association[k][i]-association[k][j]);
            }
            temp2.push_back(exp(-(par1/par)*sum));
        }
        similarity_cellline.push_back(temp2);
        temp2.clear();
    }

    for(int i=0;i<990;i++)
        for(int j=0;j<990;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<962;c++)
            if(cellline_name_962[c]==cellline_name[i])
                x=c;
        for(int d=0;d<962;d++)
            if(cellline_name_962[d]==cellline_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_cellline_962[x][y]>similarity_cellline[i][j]))
        {
            similarity_cellline[i][j]=similarity_cellline_962[x][y];
        }
    }
}

void moudle_drug(int seed_cellline)
{
    for(int i=0;i<251;i++)
        cluster_update_drug[i]=-1;
    int clunum=0;
    vector <int> sup;
    sup.clear();
    for(int i=0;i<251;i++)
    {
        if(association[i][seed_cellline]==1)
        {
            cluster_update_drug[i]=clunum;
            clunum++;
            sup.push_back(i);
        }
    }

    if(clunum!=0)
    {
    vector <int> kl;
    kl.clear();
    for(int i=0;i<251;i++)
        if(association[i][seed_cellline]==0)
        kl.push_back(i);
    int avenum=0;
    if (251%(251-kl.size())==0)
        avenum=251/(251-kl.size());
    else
        avenum=(251/(251-kl.size()))+1;
    vector <int> sup_temp;
    sup_temp.clear();
    while(1)
    {
        if(kl.size()!=0)
        {
            for(int i=0;i<kl.size();i++)
            {
            float min=36;//
            for(int j=0;j<sup.size();j++)
            {
                int flag=0;
                for(int w=0;w<sup_temp.size();w++)
                    if(j==sup_temp[w])
                    flag=1;
                if(flag==0){
                if(min>=(1-similarity_drug[kl[i]][sup[j]]))
                {
                    min=1-similarity_drug[kl[i]][sup[j]];
                    cluster_update_drug[kl[i]]=j;
                }
                }
            }
            }
            vector <int> k_l;
            k_l.clear();
            for(int j=0;j<sup.size();j++)
            {
                int flag=0;
                for(int w=0;w<sup_temp.size();w++)
                    if(j==sup_temp[w])
                    flag=1;
                if(flag==0)
                {
                vector <St> kl_temp;
                kl_temp.clear();
                for(int i=0;i<251;i++)
                {
                    if((cluster_update_drug[i]==j)&&(i!=sup[j]))
                    {
                        St s;
                        s.id=i;
                        s.sim=1-similarity_drug[i][sup[j]];
                        kl_temp.push_back(s);
                    }
                }
                if(kl_temp.size()>=avenum-1)
                {
                    for(int k=0;k<kl_temp.size()-1;k++)
                    {
                        for(int r=0;r<kl_temp.size()-k-1;r++)
                        {
                           if(kl_temp[r].sim>kl_temp[r+1].sim){
                           float tmp = kl_temp[r].sim;
                           kl_temp[r].sim = kl_temp[r+1].sim;
                           kl_temp[r+1].sim= tmp;
                           int temp_id=kl_temp[r].id;
                           kl_temp[r].id = kl_temp[r+1].id;
                           kl_temp[r+1].id= temp_id;
                        }
                        }
                    }
                    for(int g=avenum-1;g<kl_temp.size();g++)
                           k_l.push_back(kl_temp[g].id);
                    sup_temp.push_back(j);
                }
            }
            }
            kl.clear();
            for(int i=0;i<k_l.size();i++)
                kl.push_back(k_l[i]);
        }
        else
            break;
    }
    }
    int numm[clunum];
    for(int i=0;i<clunum;i++)
    {
        numm[i]=0;
        cluster c;
        c.id=i;
        for(int j=0;j<251;j++)
        {
            if(cluster_update_drug[j]==i){
                c.points.push_back(j);
                numm[i]++;
            }
        }
        clusters_drug.push_back(c);
    }
}

void moudle_cellline(int seed_drug)
{
    for(int i=0;i<990;i++)
        cluster_update_cellline[i]=-1;
    int clunum=0;
    vector <int> sup;
    sup.clear();
    for(int i=0;i<990;i++)
    {
        if(association[seed_drug][i]==1)
        {
            cluster_update_cellline[i]=clunum;
            clunum++;
            sup.push_back(i);
        }
    }
    if(clunum!=0)
    {
    vector <int> kl;
    kl.clear();
    for(int i=0;i<990;i++)
        if(association[seed_drug][i]==0)
        kl.push_back(i);
    int avenum=0;
    if (990%(990-kl.size())==0)
        avenum=990/(990-kl.size());
    else
        avenum=(990/(990-kl.size()))+1;
    vector <int> sup_temp;
    sup_temp.clear();
    while(1)
    {
        if(kl.size()!=0)
        {
            for(int i=0;i<kl.size();i++)
            {
            float min=36;//
            for(int j=0;j<sup.size();j++)
            {
                int flag=0;
                for(int w=0;w<sup_temp.size();w++)
                    if(j==sup_temp[w])
                    flag=1;
                if(flag==0){
                if(min>=(1-similarity_cellline[kl[i]][sup[j]]))
                {
                    min=1-similarity_cellline[kl[i]][sup[j]];
                    cluster_update_cellline[kl[i]]=j;
                }
                }
            }
            }
            vector <int> k_l;
            k_l.clear();
            for(int j=0;j<sup.size();j++)
            {
                int flag=0;
                for(int w=0;w<sup_temp.size();w++)
                    if(j==sup_temp[w])
                    flag=1;
                if(flag==0)
                {
                vector <St> kl_temp;
                kl_temp.clear();
                for(int i=0;i<990;i++)
                {
                    if((cluster_update_cellline[i]==j)&&(i!=sup[j]))
                    {
                        St s;
                        s.id=i;
                        s.sim=1-similarity_cellline[i][sup[j]];
                        kl_temp.push_back(s);
                    }
                }
                if(kl_temp.size()>=avenum-1)
                {
                    for(int k=0;k<kl_temp.size()-1;k++)
                    {
                        for(int r=0;r<kl_temp.size()-k-1;r++)
                        {
                           if(kl_temp[r].sim>kl_temp[r+1].sim){
                           float tmp = kl_temp[r].sim;
                           kl_temp[r].sim = kl_temp[r+1].sim;
                           kl_temp[r+1].sim= tmp;
                           int temp_id=kl_temp[r].id;
                           kl_temp[r].id = kl_temp[r+1].id;
                           kl_temp[r+1].id= temp_id;
                        }
                        }
                    }
                    for(int g=avenum-1;g<kl_temp.size();g++)
                           k_l.push_back(kl_temp[g].id);
                    sup_temp.push_back(j);
                }
            }
            }
            kl.clear();
            for(int i=0;i<k_l.size();i++)
                kl.push_back(k_l[i]);
        }
        else
            break;
    }
    }
    int numm[clunum];
    for(int i=0;i<clunum;i++)
    {
        numm[i]=0;
        cluster c;
        c.id=i;
        for(int j=0;j<990;j++)
        {
            if(cluster_update_cellline[j]==i){
                c.points.push_back(j);
                numm[i]++;
            }
        }
        clusters_cellline.push_back(c);
    }
}

float pair_score(int drug,int cellline)
{
    int degree_drug=0;
    for(int i=0;i<990;i++)
        degree_drug=degree_drug+association[drug][i];
    int degree_cellline=0;
    for(int i=0;i<251;i++)
        degree_cellline=degree_cellline+association[i][cellline];
    if((degree_drug!=0)&&(degree_cellline!=0))
    {
        moudle_drug(cellline);
        moudle_cellline(drug);
            float every[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_cellline[i].points.size();j++)
                    if((association[drug][clusters_cellline[i].points[j]]==0)&&(max<similarity_cellline[cellline][clusters_cellline[i].points[j]]))
                       max=similarity_cellline[cellline][clusters_cellline[i].points[j]];
                every[i]=max;
            }
            float every_weight[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_cellline[i].points.size();a++)
                    aver=aver+similarity_cellline[cellline][clusters_cellline[i].points[a]];
                aver=aver/clusters_cellline[i].points.size();
                every_weight[i]=aver;
            }
            float bet_cellline=0;
            float sum=0;
            float sum1=0;
            for(int i=0;i<clusters_cellline.size();i++)
            {
                sum=sum+every[i]*every_weight[i];
                sum1=sum1+every_weight[i];
            }
            bet_cellline=sum/sum1;

            float every1[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_drug[i].points.size();j++)
                    if((association[clusters_drug[i].points[j]][cellline]==0)&&(max<similarity_drug[drug][clusters_drug[i].points[j]]))
                       max=similarity_drug[drug][clusters_drug[i].points[j]];
                every1[i]=max;
            }
            float every_weight1[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_drug[i].points.size();a++)
                    aver=aver+similarity_drug[drug][clusters_drug[i].points[a]];
                aver=aver/clusters_drug[i].points.size();
                every_weight1[i]=aver;
            }
            float bet_drug=0;
            float sum2=0;
            float sum3=0;
            for(int i=0;i<clusters_drug.size();i++)
            {
                sum2=sum2+every1[i]*every_weight1[i];
                sum3=sum3+every_weight1[i];
            }
            bet_drug=sum2/sum3;

            float every2[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_cellline[i].points.size();j++)
                    if((association[drug][clusters_cellline[i].points[j]]==1)&&(max<similarity_cellline[cellline][clusters_cellline[i].points[j]]))
                       max=similarity_cellline[cellline][clusters_cellline[i].points[j]];
                every2[i]=max;
            }
            float every_weight2[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_cellline[i].points.size();a++)
                    aver=aver+similarity_cellline[cellline][clusters_cellline[i].points[a]];
                aver=aver/clusters_cellline[i].points.size();
                every_weight2[i]=aver;
            }
            float within_cellline=0;
            float sum4=0;
            float sum5=0;
            for(int i=0;i<clusters_cellline.size();i++)
            {
                sum4=sum4+every2[i]*every_weight2[i];
                sum5=sum5+every_weight2[i];
            }
            within_cellline=sum4/sum5;

            float every3[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_drug[i].points.size();j++)
                    if((association[clusters_drug[i].points[j]][cellline]==1)&&(max<similarity_drug[drug][clusters_drug[i].points[j]]))
                       max=similarity_drug[drug][clusters_drug[i].points[j]];
                every3[i]=max;
            }
            float every_weight3[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_drug[i].points.size();a++)
                    aver=aver+similarity_drug[drug][clusters_drug[i].points[a]];
                aver=aver/clusters_drug[i].points.size();
                every_weight3[i]=aver;
            }
            float within_drug=0;
            float sum6=0;
            float sum7=0;
            for(int i=0;i<clusters_drug.size();i++)
            {
                sum6=sum6+every3[i]*every_weight3[i];
                sum7=sum7+every_weight3[i];
            }
            within_drug=sum6/sum7;
            for(int s=0;s<clusters_cellline.size();s++)
                    clusters_cellline[s].points.clear();
            clusters_cellline.clear();
            for(int s=0;s<clusters_drug.size();s++)
                    clusters_drug[s].points.clear();
            clusters_drug.clear();
            return within_cellline*within_drug/(bet_cellline*bet_drug);
    }
    else if(degree_drug==0)
    {
        moudle_drug(cellline);
        float every1[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_drug[i].points.size();j++)
                    if((association[clusters_drug[i].points[j]][cellline]==1)&&(max<similarity_drug[drug][clusters_drug[i].points[j]]))
                       max=similarity_drug[drug][clusters_drug[i].points[j]];
                every1[i]=max;
            }
            float every_weight1[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_drug[i].points.size();a++)
                    aver=aver+similarity_drug[drug][clusters_drug[i].points[a]];
                aver=aver/clusters_drug[i].points.size();
                every_weight1[i]=aver;
            }
            float within_drug=0;
            float sum2=0;
            float sum3=0;
            for(int i=0;i<clusters_drug.size();i++)
            {
                sum2=sum2+every1[i]*every_weight1[i];
                sum3=sum3+every_weight1[i];
            }
            within_drug=sum2/sum3;

            float every2[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_drug[i].points.size();j++)
                    if((association[clusters_drug[i].points[j]][cellline]==0)&&(max<similarity_drug[drug][clusters_drug[i].points[j]]))
                       max=similarity_drug[drug][clusters_drug[i].points[j]];
                every2[i]=max;
            }
            float every_weight2[clusters_drug.size()];
            for(int i=0;i<clusters_drug.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_drug[i].points.size();a++)
                    aver=aver+similarity_drug[drug][clusters_drug[i].points[a]];
                aver=aver/clusters_drug[i].points.size();
                every_weight2[i]=aver;
            }
            float bet_drug=0;
            float sum4=0;
            float sum5=0;
            for(int i=0;i<clusters_drug.size();i++)
            {
                sum4=sum4+every2[i]*every_weight2[i];
                sum5=sum5+every_weight2[i];
            }
            bet_drug=sum4/sum5;

            for(int s=0;s<clusters_cellline.size();s++)
                    clusters_cellline[s].points.clear();
            clusters_cellline.clear();
            for(int s=0;s<clusters_drug.size();s++)
                    clusters_drug[s].points.clear();
            clusters_drug.clear();
            return within_drug/bet_drug;
    }
    else if(degree_cellline==0)
    {
        moudle_cellline(drug);
        float every[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_cellline[i].points.size();j++)
                    if((association[drug][clusters_cellline[i].points[j]]==1)&&(max<similarity_cellline[cellline][clusters_cellline[i].points[j]]))
                       max=similarity_cellline[cellline][clusters_cellline[i].points[j]];
                every[i]=max;
            }
            float every_weight[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_cellline[i].points.size();a++)
                    aver=aver+similarity_cellline[cellline][clusters_cellline[i].points[a]];
                aver=aver/clusters_cellline[i].points.size();
                every_weight[i]=aver;
            }
            float within_cellline=0;
            float sum=0;
            float sum1=0;
            for(int i=0;i<clusters_cellline.size();i++)
            {
                sum=sum+every[i]*every_weight[i];
                sum1=sum1+every_weight[i];
            }
            within_cellline=sum/sum1;

            float every1[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float max=0;
                for(int j=0;j<clusters_cellline[i].points.size();j++)
                    if((association[drug][clusters_cellline[i].points[j]]==0)&&(max<similarity_cellline[cellline][clusters_cellline[i].points[j]]))
                       max=similarity_cellline[cellline][clusters_cellline[i].points[j]];
                every1[i]=max;
            }
            float every_weight1[clusters_cellline.size()];
            for(int i=0;i<clusters_cellline.size();i++)
            {
                float aver=0;
                for(int a=0;a<clusters_cellline[i].points.size();a++)
                    aver=aver+similarity_cellline[cellline][clusters_cellline[i].points[a]];
                aver=aver/clusters_cellline[i].points.size();
                every_weight1[i]=aver;
            }
            float bet_cellline=0;
            float sum2=0;
            float sum3=0;
            for(int i=0;i<clusters_cellline.size();i++)
            {
                sum2=sum2+every1[i]*every_weight1[i];
                sum3=sum3+every_weight1[i];
            }
            bet_cellline=sum2/sum3;

            for(int s=0;s<clusters_cellline.size();s++)
                    clusters_cellline[s].points.clear();
            clusters_cellline.clear();
            for(int s=0;s<clusters_drug.size();s++)
                    clusters_drug[s].points.clear();
            clusters_drug.clear();
            return within_cellline/bet_cellline;
    }
}

void upload_similarity_drug(int drug)
{
    float a=0;
    for(int i=0;i<251;i++)
        for(int j=0;j<990;j++)
            a=a+association[i][j]*association[i][j];
    float par=a/251;
    for(int j=0;j<251;j++)
    {
        float sumand=0;
        for(int k=0;k<990;k++)
        {
        if((association[drug][k]==1)&&(association[j][k]==1))
            sumand++;
        }
        float par1=1/(sumand+1);
        float sum=0;
        for(int k=0;k<990;k++)
        {
            sum=sum+(association[drug][k]-association[j][k])*(association[drug][k]-association[j][k]);
        }
        similarity_drug[drug][j]=exp(-(par1/par)*sum);
        similarity_drug[j][drug]=exp(-(par1/par)*sum);
    }
    for(int j=0;j<251;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<187;c++)
            if(drug_name_187[c]==drug_name[drug])
                x=c;
        for(int d=0;d<187;d++)
            if(drug_name_187[d]==drug_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_drug_187[x][y]>similarity_drug[drug][j]))
        {
            similarity_drug[drug][j]=similarity_drug_187[x][y];
            similarity_drug[j][drug]=similarity_drug_187[x][y];
        }
    }

    for(int j=0;j<251;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<250;c++)
            if(drug_name_target[c]==drug_name[drug])
                x=c;
        for(int d=0;d<250;d++)
            if(drug_name_target[d]==drug_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_drug_target[x][y]>similarity_drug[drug][j]))
        {
            similarity_drug[drug][j]=similarity_drug_target[x][y];
            similarity_drug[j][drug]=similarity_drug_target[x][y];
        }
    }
}

void upload_similarity_cellline(int cellline)
{
    float a=0;
    for(int i=0;i<251;i++)
        for(int j=0;j<990;j++)
            a=a+association[i][j]*association[i][j];
    float par=a/990;
    for(int j=0;j<990;j++)
        {
            float sumand=0;
            for(int k=0;k<251;k++)
            {
            if((association[k][cellline]==1)&&(association[k][j]==1))
                sumand++;
            }
            float par1=1/(sumand+1);
            float sum=0;
            for(int k=0;k<251;k++)
            {
                sum=sum+(association[k][cellline]-association[k][j])*(association[k][cellline]-association[k][j]);
            }
            similarity_cellline[cellline][j]=exp(-(par1/par)*sum);
            similarity_cellline[j][cellline]=exp(-(par1/par)*sum);
        }
    for(int j=0;j<990;j++)
    {
        int x=-1;
        int y=-1;
        for(int c=0;c<962;c++)
            if(cellline_name_962[c]==cellline_name[cellline])
                x=c;
        for(int d=0;d<962;d++)
            if(cellline_name_962[d]==cellline_name[j])
                y=d;
        if((x!=-1)&&(y!=-1)&&(similarity_cellline_962[x][y]>similarity_cellline[cellline][j]))
        {
            similarity_cellline[cellline][j]=similarity_cellline_962[x][y];
            similarity_cellline[j][cellline]=similarity_cellline_962[x][y];
        }
    }
}

int main()
{
    input();
    cal_similarity_drug();
    cal_similarity_cellline();

    vector < vector <float> > score;
    ofstream out3("output data\\unloocv.txt");
    for(int i=0;i<251;i++)
    {
        vector <float> score_temp;
        for(int j=0;j<990;j++)
        {
            float ss=pair_score(i,j);
            score_temp.push_back(ss);
            out3<<ss<<" ";
        }
        score.push_back(score_temp);
        score_temp.clear();
        out3<<endl;
    }
    out3.close();

    for(int i=0;i<251;i++)
    {
        for(int j=0;j<990;j++)
        {
            if(association[i][j]==1)
            {
                association[i][j]=0;
                float temp_drug[251];
                float temp_cellline[990];
                for(int a=0;a<251;a++)
                    temp_drug[a]=similarity_drug[i][a];
                for(int a=0;a<990;a++)
                    temp_cellline[a]=similarity_cellline[j][a];
                upload_similarity_drug(i);
                upload_similarity_cellline(j);
                score[i][j]=pair_score(i,j);
                association[i][j]=1;
                for(int a=0;a<251;a++)
                {
                    similarity_drug[i][a]=temp_drug[a];
                    similarity_drug[a][i]=temp_drug[a];
                }
                for(int a=0;a<990;a++)
                {
                    similarity_cellline[j][a]=temp_cellline[a];
                    similarity_cellline[a][j]=temp_cellline[a];
                }
            }
        }
    }
    ofstream out4("output data\\loocv.txt");
    for(int i=0;i<251;i++)
    {
        for(int j=0;j<990;j++)
            out4<<score[i][j]<<" ";
        out4<<endl;
    }
    out4.close();
    drug_name_target.clear();
    for(int a=0;a<251;a++)
        association[a].clear();
    association.clear();
    for(int a=0;a<251;a++)
        score[a].clear();
    score.clear();
    for(int a=0;a<962;a++)
        similarity_cellline_962[a].clear();
    similarity_cellline_962.clear();
    for(int a=0;a<251;a++)
        similarity_drug[a].clear();
    similarity_drug.clear();
    for(int a=0;a<250;a++)
        similarity_drug_target[a].clear();
    similarity_drug_target.clear();
    for(int a=0;a<990;a++)
        similarity_cellline[a].clear();
    similarity_cellline.clear();
    return 0;
}
