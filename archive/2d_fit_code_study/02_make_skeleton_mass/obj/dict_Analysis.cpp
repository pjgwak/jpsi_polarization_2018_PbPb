// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME objdIdict_Analysis
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
#include "src/ModelBuilder.h"
#include "src/JpsiFitter.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *ModelBuilder_Dictionary();
   static void ModelBuilder_TClassManip(TClass*);
   static void delete_ModelBuilder(void *p);
   static void deleteArray_ModelBuilder(void *p);
   static void destruct_ModelBuilder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModelBuilder*)
   {
      ::ModelBuilder *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ModelBuilder));
      static ::ROOT::TGenericClassInfo 
         instance("ModelBuilder", "src/ModelBuilder.h", 8,
                  typeid(::ModelBuilder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ModelBuilder_Dictionary, isa_proxy, 4,
                  sizeof(::ModelBuilder) );
      instance.SetDelete(&delete_ModelBuilder);
      instance.SetDeleteArray(&deleteArray_ModelBuilder);
      instance.SetDestructor(&destruct_ModelBuilder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModelBuilder*)
   {
      return GenerateInitInstanceLocal(static_cast<::ModelBuilder*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ModelBuilder*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ModelBuilder_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ModelBuilder*>(nullptr))->GetClass();
      ModelBuilder_TClassManip(theClass);
   return theClass;
   }

   static void ModelBuilder_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *JpsiFitter_Dictionary();
   static void JpsiFitter_TClassManip(TClass*);
   static void delete_JpsiFitter(void *p);
   static void deleteArray_JpsiFitter(void *p);
   static void destruct_JpsiFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JpsiFitter*)
   {
      ::JpsiFitter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JpsiFitter));
      static ::ROOT::TGenericClassInfo 
         instance("JpsiFitter", "src/JpsiFitter.h", 7,
                  typeid(::JpsiFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JpsiFitter_Dictionary, isa_proxy, 4,
                  sizeof(::JpsiFitter) );
      instance.SetDelete(&delete_JpsiFitter);
      instance.SetDeleteArray(&deleteArray_JpsiFitter);
      instance.SetDestructor(&destruct_JpsiFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JpsiFitter*)
   {
      return GenerateInitInstanceLocal(static_cast<::JpsiFitter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::JpsiFitter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JpsiFitter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::JpsiFitter*>(nullptr))->GetClass();
      JpsiFitter_TClassManip(theClass);
   return theClass;
   }

   static void JpsiFitter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ModelBuilder(void *p) {
      delete (static_cast<::ModelBuilder*>(p));
   }
   static void deleteArray_ModelBuilder(void *p) {
      delete [] (static_cast<::ModelBuilder*>(p));
   }
   static void destruct_ModelBuilder(void *p) {
      typedef ::ModelBuilder current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ModelBuilder

namespace ROOT {
   // Wrapper around operator delete
   static void delete_JpsiFitter(void *p) {
      delete (static_cast<::JpsiFitter*>(p));
   }
   static void deleteArray_JpsiFitter(void *p) {
      delete [] (static_cast<::JpsiFitter*>(p));
   }
   static void destruct_JpsiFitter(void *p) {
      typedef ::JpsiFitter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::JpsiFitter

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 428,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *maplEstringcOvectorlEdoublegRsPgR_Dictionary();
   static void maplEstringcOvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlEdoublegRsPgR(void *p = nullptr);
   static void *newArray_maplEstringcOvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlEdoublegRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlEdoublegRsPgR(void *p);
   static void destruct_maplEstringcOvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<double> >*)
   {
      map<string,vector<double> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<double> >", -2, "map", 102,
                  typeid(map<string,vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,vector<double> >) );
      instance.SetNew(&new_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<double> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,vector<double> >","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<double, std::allocator<double> >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::vector<double, std::allocator<double> > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,vector<double> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,vector<double> >*>(nullptr))->GetClass();
      maplEstringcOvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,vector<double> > : new map<string,vector<double> >;
   }
   static void *newArray_maplEstringcOvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,vector<double> >[nElements] : new map<string,vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      delete (static_cast<map<string,vector<double> >*>(p));
   }
   static void deleteArray_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      delete [] (static_cast<map<string,vector<double> >*>(p));
   }
   static void destruct_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      typedef map<string,vector<double> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<double> >

namespace {
  void TriggerDictionaryInitialization_dict_Analysis_Impl() {
    static const char* headers[] = {
"src/ModelBuilder.h",
"src/JpsiFitter.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/hep319/anaconda3/envs/root632/include/",
"/work/pjgwak/pol24/archive/2d_fit_code_study/02_make_skeleton_mass/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict_Analysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/ModelBuilder.h")))  ModelBuilder;
class __attribute__((annotate("$clingAutoload$src/JpsiFitter.h")))  JpsiFitter;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict_Analysis dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/ModelBuilder.h"
#include "src/JpsiFitter.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"JpsiFitter", payloadCode, "@",
"ModelBuilder", payloadCode, "@",
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
