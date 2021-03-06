#ifndef CONFIG_H
#define CONFIG_H


#include <vector>
#include <string>
#include <utility>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <Handler.h>

//#define xmlFreeDoc delete


using namespace std ;
 

class Configurator;

class Configurable{ // base class - abstract
 public:
  virtual void Clear()=0;
  virtual void Init()=0;
  virtual void Config(Configurator &)=0;
  static string getElementContent (xmlDocPtr doc, const char * key, const xmlNode * node);
  static string getElementContent (Configurator&c, const char * key, const xmlNode * node);
  static vector<string> getElementVector (Configurator&c, const char * key, const xmlNode * node);
  static vector<string> getElementVector (xmlDocPtr doc, const char * key, const xmlNode * node);
  static vector<pair<string, string> > getNodeContentList (xmlDocPtr doc, const xmlNode * node) ;
  static vector<pair<string, string> > getNodeContentList (Configurator & c, const xmlNode * node) ;


};

class Configurator 
{
 private:
  //---
 public:
  xmlDoc *doc;// = NULL;
  xmlNode *root_element;// = NULL;
  string xmlFileName;
  // --- constructor
  Configurator();
  // --- destructor
  ~Configurator();
  // --- Initialization 
  virtual void Init();
  virtual void Clear();
  void ReadFromStream(string stream);
  static inline double GetDouble(string value) { return atof(value.c_str());}
  static int    GetInt   (string value);// { return atoi(value.c_str());}
  static inline long   GetLong  (string value) { return atol(value.c_str());}
  static void  GetVecInt(string value,vector<int> &v) ;
  static void  GetVecFloat(string value,vector<float> &v) ;
  
  // -------------XML
};

#endif
