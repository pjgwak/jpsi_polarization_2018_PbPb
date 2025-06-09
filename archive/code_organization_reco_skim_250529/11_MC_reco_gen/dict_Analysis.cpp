// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dict_Analysis
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "EventLoop.h"
#include "Data.h"
#include "Algorithm.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *Data_Dictionary();
   static void Data_TClassManip(TClass*);
   static void delete_Data(void *p);
   static void deleteArray_Data(void *p);
   static void destruct_Data(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Data*)
   {
      ::Data *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Data));
      static ::ROOT::TGenericClassInfo 
         instance("Data", "Data.h", 7,
                  typeid(::Data), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Data_Dictionary, isa_proxy, 4,
                  sizeof(::Data) );
      instance.SetDelete(&delete_Data);
      instance.SetDeleteArray(&deleteArray_Data);
      instance.SetDestructor(&destruct_Data);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Data*)
   {
      return GenerateInitInstanceLocal(static_cast<::Data*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Data*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Data_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Data*>(nullptr))->GetClass();
      Data_TClassManip(theClass);
   return theClass;
   }

   static void Data_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DataRun2_Dictionary();
   static void DataRun2_TClassManip(TClass*);
   static void delete_DataRun2(void *p);
   static void deleteArray_DataRun2(void *p);
   static void destruct_DataRun2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DataRun2*)
   {
      ::DataRun2 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DataRun2));
      static ::ROOT::TGenericClassInfo 
         instance("DataRun2", "Data.h", 59,
                  typeid(::DataRun2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DataRun2_Dictionary, isa_proxy, 4,
                  sizeof(::DataRun2) );
      instance.SetDelete(&delete_DataRun2);
      instance.SetDeleteArray(&deleteArray_DataRun2);
      instance.SetDestructor(&destruct_DataRun2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DataRun2*)
   {
      return GenerateInitInstanceLocal(static_cast<::DataRun2*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::DataRun2*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DataRun2_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::DataRun2*>(nullptr))->GetClass();
      DataRun2_TClassManip(theClass);
   return theClass;
   }

   static void DataRun2_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DataRun3_Dictionary();
   static void DataRun3_TClassManip(TClass*);
   static void delete_DataRun3(void *p);
   static void deleteArray_DataRun3(void *p);
   static void destruct_DataRun3(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DataRun3*)
   {
      ::DataRun3 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DataRun3));
      static ::ROOT::TGenericClassInfo 
         instance("DataRun3", "Data.h", 77,
                  typeid(::DataRun3), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DataRun3_Dictionary, isa_proxy, 4,
                  sizeof(::DataRun3) );
      instance.SetDelete(&delete_DataRun3);
      instance.SetDeleteArray(&deleteArray_DataRun3);
      instance.SetDestructor(&destruct_DataRun3);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DataRun3*)
   {
      return GenerateInitInstanceLocal(static_cast<::DataRun3*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::DataRun3*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DataRun3_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::DataRun3*>(nullptr))->GetClass();
      DataRun3_TClassManip(theClass);
   return theClass;
   }

   static void DataRun3_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Algorithm_Dictionary();
   static void Algorithm_TClassManip(TClass*);
   static void *new_Algorithm(void *p = nullptr);
   static void *newArray_Algorithm(Long_t size, void *p);
   static void delete_Algorithm(void *p);
   static void deleteArray_Algorithm(void *p);
   static void destruct_Algorithm(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Algorithm*)
   {
      ::Algorithm *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Algorithm));
      static ::ROOT::TGenericClassInfo 
         instance("Algorithm", "Algorithm.h", 10,
                  typeid(::Algorithm), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Algorithm_Dictionary, isa_proxy, 4,
                  sizeof(::Algorithm) );
      instance.SetNew(&new_Algorithm);
      instance.SetNewArray(&newArray_Algorithm);
      instance.SetDelete(&delete_Algorithm);
      instance.SetDeleteArray(&deleteArray_Algorithm);
      instance.SetDestructor(&destruct_Algorithm);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Algorithm*)
   {
      return GenerateInitInstanceLocal(static_cast<::Algorithm*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Algorithm*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Algorithm_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Algorithm*>(nullptr))->GetClass();
      Algorithm_TClassManip(theClass);
   return theClass;
   }

   static void Algorithm_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *EventLoop_Dictionary();
   static void EventLoop_TClassManip(TClass*);
   static void *new_EventLoop(void *p = nullptr);
   static void *newArray_EventLoop(Long_t size, void *p);
   static void delete_EventLoop(void *p);
   static void deleteArray_EventLoop(void *p);
   static void destruct_EventLoop(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventLoop*)
   {
      ::EventLoop *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EventLoop));
      static ::ROOT::TGenericClassInfo 
         instance("EventLoop", "EventLoop.h", 9,
                  typeid(::EventLoop), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &EventLoop_Dictionary, isa_proxy, 4,
                  sizeof(::EventLoop) );
      instance.SetNew(&new_EventLoop);
      instance.SetNewArray(&newArray_EventLoop);
      instance.SetDelete(&delete_EventLoop);
      instance.SetDeleteArray(&deleteArray_EventLoop);
      instance.SetDestructor(&destruct_EventLoop);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventLoop*)
   {
      return GenerateInitInstanceLocal(static_cast<::EventLoop*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::EventLoop*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *EventLoop_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::EventLoop*>(nullptr))->GetClass();
      EventLoop_TClassManip(theClass);
   return theClass;
   }

   static void EventLoop_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Data(void *p) {
      delete (static_cast<::Data*>(p));
   }
   static void deleteArray_Data(void *p) {
      delete [] (static_cast<::Data*>(p));
   }
   static void destruct_Data(void *p) {
      typedef ::Data current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Data

namespace ROOT {
   // Wrapper around operator delete
   static void delete_DataRun2(void *p) {
      delete (static_cast<::DataRun2*>(p));
   }
   static void deleteArray_DataRun2(void *p) {
      delete [] (static_cast<::DataRun2*>(p));
   }
   static void destruct_DataRun2(void *p) {
      typedef ::DataRun2 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::DataRun2

namespace ROOT {
   // Wrapper around operator delete
   static void delete_DataRun3(void *p) {
      delete (static_cast<::DataRun3*>(p));
   }
   static void deleteArray_DataRun3(void *p) {
      delete [] (static_cast<::DataRun3*>(p));
   }
   static void destruct_DataRun3(void *p) {
      typedef ::DataRun3 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::DataRun3

namespace ROOT {
   // Wrappers around operator new
   static void *new_Algorithm(void *p) {
      return  p ? new(p) ::Algorithm : new ::Algorithm;
   }
   static void *newArray_Algorithm(Long_t nElements, void *p) {
      return p ? new(p) ::Algorithm[nElements] : new ::Algorithm[nElements];
   }
   // Wrapper around operator delete
   static void delete_Algorithm(void *p) {
      delete (static_cast<::Algorithm*>(p));
   }
   static void deleteArray_Algorithm(void *p) {
      delete [] (static_cast<::Algorithm*>(p));
   }
   static void destruct_Algorithm(void *p) {
      typedef ::Algorithm current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Algorithm

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventLoop(void *p) {
      return  p ? new(p) ::EventLoop : new ::EventLoop;
   }
   static void *newArray_EventLoop(Long_t nElements, void *p) {
      return p ? new(p) ::EventLoop[nElements] : new ::EventLoop[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventLoop(void *p) {
      delete (static_cast<::EventLoop*>(p));
   }
   static void deleteArray_EventLoop(void *p) {
      delete [] (static_cast<::EventLoop*>(p));
   }
   static void destruct_EventLoop(void *p) {
      typedef ::EventLoop current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::EventLoop

namespace ROOT {
   static TClass *vectorlETStringgR_Dictionary();
   static void vectorlETStringgR_TClassManip(TClass*);
   static void *new_vectorlETStringgR(void *p = nullptr);
   static void *newArray_vectorlETStringgR(Long_t size, void *p);
   static void delete_vectorlETStringgR(void *p);
   static void deleteArray_vectorlETStringgR(void *p);
   static void destruct_vectorlETStringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TString>*)
   {
      vector<TString> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TString>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TString>", -2, "vector", 428,
                  typeid(vector<TString>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETStringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TString>) );
      instance.SetNew(&new_vectorlETStringgR);
      instance.SetNewArray(&newArray_vectorlETStringgR);
      instance.SetDelete(&delete_vectorlETStringgR);
      instance.SetDeleteArray(&deleteArray_vectorlETStringgR);
      instance.SetDestructor(&destruct_vectorlETStringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TString> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TString>","std::vector<TString, std::allocator<TString> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TString>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETStringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TString>*>(nullptr))->GetClass();
      vectorlETStringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETStringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETStringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TString> : new vector<TString>;
   }
   static void *newArray_vectorlETStringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TString>[nElements] : new vector<TString>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETStringgR(void *p) {
      delete (static_cast<vector<TString>*>(p));
   }
   static void deleteArray_vectorlETStringgR(void *p) {
      delete [] (static_cast<vector<TString>*>(p));
   }
   static void destruct_vectorlETStringgR(void *p) {
      typedef vector<TString> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TString>

namespace ROOT {
   static TClass *vectorlEAlgorithmmUgR_Dictionary();
   static void vectorlEAlgorithmmUgR_TClassManip(TClass*);
   static void *new_vectorlEAlgorithmmUgR(void *p = nullptr);
   static void *newArray_vectorlEAlgorithmmUgR(Long_t size, void *p);
   static void delete_vectorlEAlgorithmmUgR(void *p);
   static void deleteArray_vectorlEAlgorithmmUgR(void *p);
   static void destruct_vectorlEAlgorithmmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Algorithm*>*)
   {
      vector<Algorithm*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Algorithm*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Algorithm*>", -2, "vector", 428,
                  typeid(vector<Algorithm*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEAlgorithmmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Algorithm*>) );
      instance.SetNew(&new_vectorlEAlgorithmmUgR);
      instance.SetNewArray(&newArray_vectorlEAlgorithmmUgR);
      instance.SetDelete(&delete_vectorlEAlgorithmmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEAlgorithmmUgR);
      instance.SetDestructor(&destruct_vectorlEAlgorithmmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Algorithm*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Algorithm*>","std::vector<Algorithm*, std::allocator<Algorithm*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Algorithm*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEAlgorithmmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Algorithm*>*>(nullptr))->GetClass();
      vectorlEAlgorithmmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEAlgorithmmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEAlgorithmmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Algorithm*> : new vector<Algorithm*>;
   }
   static void *newArray_vectorlEAlgorithmmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Algorithm*>[nElements] : new vector<Algorithm*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEAlgorithmmUgR(void *p) {
      delete (static_cast<vector<Algorithm*>*>(p));
   }
   static void deleteArray_vectorlEAlgorithmmUgR(void *p) {
      delete [] (static_cast<vector<Algorithm*>*>(p));
   }
   static void destruct_vectorlEAlgorithmmUgR(void *p) {
      typedef vector<Algorithm*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Algorithm*>

namespace {
  void TriggerDictionaryInitialization_dict_Analysis_Impl() {
    static const char* headers[] = {
"EventLoop.h",
"Data.h",
"Algorithm.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/hep319/anaconda3/envs/root632/include/",
"/work/pjgwak/pol24/archive/code_organization_reco_skim_250529/10_run2_data_added/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict_Analysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Data.h")))  __attribute__((annotate("$clingAutoload$EventLoop.h")))  Data;
class __attribute__((annotate("$clingAutoload$Data.h")))  __attribute__((annotate("$clingAutoload$EventLoop.h")))  DataRun2;
class __attribute__((annotate("$clingAutoload$Data.h")))  __attribute__((annotate("$clingAutoload$EventLoop.h")))  DataRun3;
class __attribute__((annotate("$clingAutoload$Algorithm.h")))  __attribute__((annotate("$clingAutoload$EventLoop.h")))  Algorithm;
class __attribute__((annotate("$clingAutoload$EventLoop.h")))  EventLoop;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict_Analysis dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "EventLoop.h"
#include "Data.h"
#include "Algorithm.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Algorithm", payloadCode, "@",
"Data", payloadCode, "@",
"DataRun2", payloadCode, "@",
"DataRun3", payloadCode, "@",
"EventLoop", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict_Analysis",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Analysis_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Analysis_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict_Analysis() {
  TriggerDictionaryInitialization_dict_Analysis_Impl();
}
