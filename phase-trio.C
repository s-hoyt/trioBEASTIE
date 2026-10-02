/****************************************************************
 phase-trio.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Essex.H"
#include "BOOM/VcfReader.H"
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 This program performs trio phasing.  All combinations can be
 perfectly phased except triple hets.  There are two possible
 ways to phase a triple het: 10 01 10, or 01 10 01 (M/F/C).  This
 program arbitrarily phases one of these two ways.  Downstream
 models should sum over these two phases by swapping the alleles
 in all three individuals to get the other phase.
 ****************************************************************/

enum Individual { MOTHER=0, FATHER=1, CHILD=2 };
enum Allele { REF=0, ALT=1 };
enum MaternalPaternal { MAT=0, PAT=1 };

class Application {
  Map<String,String> phasingMap;
  void initMap(Map<String,String> &);
  bool phase(Genotype &mother,Genotype &father,Genotype &child);
  Genotype getEssexGT(Essex::Node *siteGenotype,String label);
  void installGT(Essex::Node *,const String &label,const Genotype &);
  int getEssexNumericChild(Essex::CompositeNode *,int whichChild);
  void setEssexNumericChild(Essex::CompositeNode *,int which,int value);
  String compactString(const Genotype &);
  void install(const String &encoded,int first,int second,Genotype &);
  void installSuccess(bool phased,Essex::Node *);
  void phaseCounts(Genotype mother,Genotype father,Genotype child,   
		   Essex::Node *site);
  void phaseCounts(Genotype G,String label,Essex::Node *counts);
  void swapCounts(Essex::Node *counts,String label);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  // ctor

  initMap(phasingMap);
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("phase-trio <in.essex> <out.essex>");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);

  // Create output file
  ofstream os(outfile);
  
  // Process each gene in the input file
  Essex::Parser parser(infile);
  Essex::Node *root;
  while(root=parser.nextElem()) {
    Vector<Essex::Node*> sites;
    //root->findDescendents("genotypes",sites);
    root->findDescendents("site",sites);
    for(Vector<Essex::Node*>::iterator cur=sites.begin(), end=sites.end() ;
	cur!=end ; ++cur) {
      Essex::Node *site=*cur;
      Essex::Node *genotypes=site->findChild("genotypes");
      Genotype mother=getEssexGT(genotypes,"mother");
      Genotype father=getEssexGT(genotypes,"father");
      Genotype child=getEssexGT(genotypes,"child");
      bool success=phase(mother,father,child);
      installGT(genotypes,"child",child);
      installGT(genotypes,"mother",mother);
      installGT(genotypes,"father",father);
      installSuccess(success,site);
      //if(success) phaseCounts(mother,father,child,site); // wrong
      phaseCounts(mother,father,child,site); // right, but lower accuracy?!
    }
    root->printOn(os); os<<endl;
  }

  return 0;
}



void Application::phaseCounts(Genotype mother,Genotype father,Genotype child,
			      Essex::Node *site)
{
  // PRECONDITION: this site was able to be phased

  Essex::Node *counts=site->findChild("counts");
  phaseCounts(mother,"mother",counts);
  phaseCounts(father,"father",counts);
  phaseCounts(child,"child",counts);
}



void Application::phaseCounts(Genotype G,String label,Essex::Node *counts)
{
  // PRECONDITION: genotype has two alleles
  // Prior to this function, the counts are ordered REF/ALT, but after
  // this function, they will be ordered MAT/PAT

  if(!G.isHet()) return;
  if(G[MAT]==ALT) swapCounts(counts,label);
}



void Application::swapCounts(Essex::Node *counts,String label)
{
  Essex::CompositeNode *node=
    dynamic_cast<Essex::CompositeNode*>(counts->findChild(label));
  if(!node) throw "Cannot cast in swapCounts";
  if(node->getNumChildren()!=2) 
    throw "Wrong number of children in swapCounts()";
  Essex::NumericNode *mat=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(0));
  Essex::NumericNode *pat=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(1));
  float temp=mat->getValue();
  mat->setValue(pat->getValue());
  pat->setValue(temp);
}



void Application::installSuccess(bool phased,Essex::Node *site)
{
  Essex::CompositeNode *comp=static_cast<Essex::CompositeNode*>(site);
  comp->append("phased",phased ? 1 : 0);
}


int Application::getEssexNumericChild(Essex::CompositeNode *node,int which)
{
  Essex::NumericNode *child=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(which));
  if(!child) throw "Child is not numeric in getEssexNumericChild";
  return child->getValue();
}



void Application::setEssexNumericChild(Essex::CompositeNode *node,int which,
				      int value)
{
  Essex::NumericNode *child=
    dynamic_cast<Essex::NumericNode*>(node->getIthChild(which));
  if(!child) throw "Child is not numeric in getEssexNumericChild";
  child->setValue(value);
}



void Application::installGT(Essex::Node *site,const String &label,
			    const Genotype &G)
{
  Essex::CompositeNode *child=
    dynamic_cast<Essex::CompositeNode*>(site->findChild(label));
  if(!child) throw label+" not found in getEssexGT";
  if(child->getNumChildren()!=2) throw label+" has wrong number of children";
  setEssexNumericChild(child,0,G[0]);
  setEssexNumericChild(child,1,G[1]);
}



Genotype Application::getEssexGT(Essex::Node *site,String label)
{
  Essex::CompositeNode *child=
    dynamic_cast<Essex::CompositeNode*>(site->findChild(label));
  if(!child) throw label+" not found in getEssexGT";
  if(child->getNumChildren()!=2) throw label+" has wrong number of children";
  Genotype g;
  int allele1=getEssexNumericChild(child,0);
  int allele2=getEssexNumericChild(child,1);
  g.addAllele(allele1);
  g.addAllele(allele2);
  return g;
}



String Application::compactString(const Genotype &G)
{
  String s="";
  const int numAlleles=G.numAlleles();
  for(int i=0 ; i<numAlleles ; ++i) s+=String(G[i]);
  if(s=="10") s="01"; // Normalize unphased genotypes for easy lookup
  return s;
}



void Application::initMap(Map<String,String> &M)
{
  // PRECONDITION: each individual has been normalized, so that 10
  // doesn't occur in the inputs (for an individual)

  // The first two columns are mother, second pair is father, third is child
  
  M["000000"]="000000";
  M["000100"]="000100";
  M["000101"]="001001";
  M["001101"]="001101";
  M["010000"]="010000";
  M["010001"]="100010";
  M["010100"]="010100";
  M["010111"]="101011";
  M["011101"]="011101";
  M["011111"]="101111";
  M["110001"]="110010";
  M["110101"]="110110";
  M["110111"]="111011";
  M["111111"]="111111";
  M["010101"]="100110"; // triple het => phase arbitrarily & marginalize later

  //M["110000"]=""; DE NOVO
  //M["000001"]=""; DE NOVO, CAN'T PHASE
  //M["000010"]=""; DE NOVO, CAN'T PHASE
  //M["000011"]=""; TWO DE NOVOS
  //M["000110"]=""; UNNORMALIZED
  //M["000111"]=""; DE NOVO
  //M["001100"]=""; DE NOVO
  //M["001000"]=""; UNNORMALIZED
  //M["001001"]=""; UNNORMALIZED
  //M["001010"]=""; UNNORMALIZED
  //M["001011"]=""; UNNORMALIZED
  //M["001110"]=""; UNNORMALIZED
  //M["001111"]=""; DE NOVO
  //M["010010"]=""; UNNORMALIZED
  //M["010011"]=""; DE NOVO
  //M["010110"]=""; UNNORMALIZED
  //M["011000"]=""; UNNORMALIZED
  //M["011001"]=""; UNNORMALIZED
  //M["011010"]=""; UNNORMALIZED
  //M["011011"]=""; UNNORMALIZED
  //M["011100"]=""; DE NOVO
  //M["011110"]=""; UNNORMALIZED
  //M["100000"]=""; UNNORMALIZED
  //M["100001"]=""; UNNORMALIZED
  //M["100010"]=""; UNNORMALIZED
  //M["100011"]=""; UNNORMALIZED
  //M["100100"]=""; UNNORMALIZED
  //M["100101"]=""; UNNORMALIZED
  //M["100110"]=""; UNNORMALIZED
  //M["100111"]=""; UNNORMALIZED
  //M["101000"]=""; UNNORMALIZED
  //M["101001"]=""; UNNORMALIZED
  //M["101010"]=""; UNNORMALIZED
  //M["101011"]=""; UNNORMALIZED
  //M["101100"]=""; UNNORMALIZED
  //M["101101"]=""; UNNORMALIZED
  //M["101110"]=""; UNNORMALIZED
  //M["101111"]=""; UNNORMALIZED
  //M["110010"]=""; UNNORMALIZED
  //M["110011"]=""; DE NOVO
  //M["110100"]=""; DE NOVO
  //M["110110"]=""; UNNORMALIZED
  //M["111000"]=""; UNNORMALIZED
  //M["111001"]=""; UNNORMALIZED
  //M["111010"]=""; UNNORMALIZED
  //M["111011"]=""; UNNORMALIZED
  //M["111100"]=""; DE NOVO
  //M["111101"]=""; DE NOVO
  //M["111110"]=""; UNNORMALIZED
}



bool Application::phase(Genotype &mother,Genotype &father,Genotype &child)
{
  // PRECONDITION: all three genotypes have two alleles each

  const String encoded=
    compactString(mother)+compactString(father)+compactString(child);
  bool canPhase=true;
  if(encoded=="010101") canPhase=false; // triple het
  if(!phasingMap.isDefined(encoded))
    throw String("Genotype encoding is not defined: ")+encoded;
  String phased=phasingMap[encoded];
  install(phased,0,1,mother);
  install(phased,2,3,father);
  install(phased,4,5,child);
  return canPhase;
}



void Application::install(const String &encoded,int first,int second,
			  Genotype &G)
{
  Vector<int> &v=G.asVector();
  const char *p=encoded.c_str();
  v[0]=p[first]-'0';
  v[1]=p[second]-'0';
}



