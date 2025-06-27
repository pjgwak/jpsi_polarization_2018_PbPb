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
#include "src/McMassFitter.h"
#include "src/MassFitter.h"
#include "src/SPlotter.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *McMassFitter_Dictionary();
   static void McMassFitter_TClassManip(TClass*);
   static void delete_McMassFitter(void *p);
   static void deleteArray_McMassFitter(void *p);
   static void destruct_McMassFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::McMassFitter*)
   {
      ::McMassFitter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::McMassFitter));
      static ::ROOT::TGenericClassInfo 
         instance("McMassFitter", "McMassFitter.h", 12,
                  typeid(::McMassFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &McMassFitter_Dictionary, isa_proxy, 4,
                  sizeof(::McMassFitter) );
      instance.SetDelete(&delete_McMassFitter);
      instance.SetDeleteArray(&deleteArray_McMassFitter);
      instance.SetDestructor(&destruct_McMassFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::McMassFitter*)
   {
      return GenerateInitInstanceLocal(static_cast<::McMassFitter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::McMassFitter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *McMassFitter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::McMassFitter*>(nullptr))->GetClass();
      McMassFitter_TClassManip(theClass);
   return theClass;
   }

   static void McMassFitter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *MassFitter_Dictionary();
   static void MassFitter_TClassManip(TClass*);
   static void delete_MassFitter(void *p);
   static void deleteArray_MassFitter(void *p);
   static void destruct_MassFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MassFitter*)
   {
      ::MassFitter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MassFitter));
      static ::ROOT::TGenericClassInfo 
         instance("MassFitter", "MassFitter.h", 13,
                  typeid(::MassFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MassFitter_Dictionary, isa_proxy, 4,
                  sizeof(::MassFitter) );
      instance.SetDelete(&delete_MassFitter);
      instance.SetDeleteArray(&deleteArray_MassFitter);
      instance.SetDestructor(&destruct_MassFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MassFitter*)
   {
      return GenerateInitInstanceLocal(static_cast<::MassFitter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MassFitter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MassFitter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::MassFitter*>(nullptr))->GetClass();
      MassFitter_TClassManip(theClass);
   return theClass;
   }

   static void MassFitter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SPlotter_Dictionary();
   static void SPlotter_TClassManip(TClass*);
   static void delete_SPlotter(void *p);
   static void deleteArray_SPlotter(void *p);
   static void destruct_SPlotter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SPlotter*)
   {
      ::SPlotter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SPlotter));
      static ::ROOT::TGenericClassInfo 
         instance("SPlotter", "SPlotter.h", 20,
                  typeid(::SPlotter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SPlotter_Dictionary, isa_proxy, 4,
                  sizeof(::SPlotter) );
      instance.SetDelete(&delete_SPlotter);
      instance.SetDeleteArray(&deleteArray_SPlotter);
      instance.SetDestructor(&destruct_SPlotter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SPlotter*)
   {
      return GenerateInitInstanceLocal(static_cast<::SPlotter*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::SPlotter*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SPlotter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::SPlotter*>(nullptr))->GetClass();
      SPlotter_TClassManip(theClass);
   return theClass;
   }

   static void SPlotter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_McMassFitter(void *p) {
      delete (static_cast<::McMassFitter*>(p));
   }
   static void deleteArray_McMassFitter(void *p) {
      delete [] (static_cast<::McMassFitter*>(p));
   }
   static void destruct_McMassFitter(void *p) {
      typedef ::McMassFitter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::McMassFitter

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MassFitter(void *p) {
      delete (static_cast<::MassFitter*>(p));
   }
   static void deleteArray_MassFitter(void *p) {
      delete [] (static_cast<::MassFitter*>(p));
   }
   static void destruct_MassFitter(void *p) {
      typedef ::MassFitter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::MassFitter

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SPlotter(void *p) {
      delete (static_cast<::SPlotter*>(p));
   }
   static void deleteArray_SPlotter(void *p) {
      delete [] (static_cast<::SPlotter*>(p));
   }
   static void destruct_SPlotter(void *p) {
      typedef ::SPlotter current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::SPlotter

namespace {
  void TriggerDictionaryInitialization_dict_Analysis_Impl() {
    static const char* headers[] = {
"src/McMassFitter.h",
"src/MassFitter.h",
"src/SPlotter.h",
nullptr
    };
    static const char* includePaths[] = {
"src",
"/home/hep319/anaconda3/envs/root632/include/",
"/work/pjgwak/pol24/archive/2d_fit_code_study/05_rework_from_mc_mass/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict_Analysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/McMassFitter.h")))  McMassFitter;
class __attribute__((annotate("$clingAutoload$src/MassFitter.h")))  MassFitter;
class __attribute__((annotate("$clingAutoload$src/SPlotter.h")))  SPlotter;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict_Analysis dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/McMassFitter.h"
#include "src/MassFitter.h"
#include "src/SPlotter.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MassFitter", payloadCode, "@",
"McMassFitter", payloadCode, "@",
"SPlotter", payloadCode, "@",
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
