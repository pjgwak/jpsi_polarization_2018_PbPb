// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME objdIdict_FlowSkim
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
#include "src/FlowSkimRun3DataPbPb.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *FlowSkimRun3DataPbPb_Dictionary();
   static void FlowSkimRun3DataPbPb_TClassManip(TClass*);
   static void delete_FlowSkimRun3DataPbPb(void *p);
   static void deleteArray_FlowSkimRun3DataPbPb(void *p);
   static void destruct_FlowSkimRun3DataPbPb(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FlowSkimRun3DataPbPb*)
   {
      ::FlowSkimRun3DataPbPb *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FlowSkimRun3DataPbPb));
      static ::ROOT::TGenericClassInfo 
         instance("FlowSkimRun3DataPbPb", "FlowSkimRun3DataPbPb.h", 12,
                  typeid(::FlowSkimRun3DataPbPb), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FlowSkimRun3DataPbPb_Dictionary, isa_proxy, 4,
                  sizeof(::FlowSkimRun3DataPbPb) );
      instance.SetDelete(&delete_FlowSkimRun3DataPbPb);
      instance.SetDeleteArray(&deleteArray_FlowSkimRun3DataPbPb);
      instance.SetDestructor(&destruct_FlowSkimRun3DataPbPb);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FlowSkimRun3DataPbPb*)
   {
      return GenerateInitInstanceLocal(static_cast<::FlowSkimRun3DataPbPb*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FlowSkimRun3DataPbPb*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FlowSkimRun3DataPbPb_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::FlowSkimRun3DataPbPb*>(nullptr))->GetClass();
      FlowSkimRun3DataPbPb_TClassManip(theClass);
   return theClass;
   }

   static void FlowSkimRun3DataPbPb_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_FlowSkimRun3DataPbPb(void *p) {
      delete (static_cast<::FlowSkimRun3DataPbPb*>(p));
   }
   static void deleteArray_FlowSkimRun3DataPbPb(void *p) {
      delete [] (static_cast<::FlowSkimRun3DataPbPb*>(p));
   }
   static void destruct_FlowSkimRun3DataPbPb(void *p) {
      typedef ::FlowSkimRun3DataPbPb current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FlowSkimRun3DataPbPb

namespace {
  void TriggerDictionaryInitialization_dict_FlowSkim_Impl() {
    static const char* headers[] = {
"src/FlowSkimRun3DataPbPb.h",
nullptr
    };
    static const char* includePaths[] = {
"src",
"/home/hep319/anaconda3/envs/root632/include/",
"/work/pjgwak/pol24/input_flowskim/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict_FlowSkim dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/FlowSkimRun3DataPbPb.h")))  FlowSkimRun3DataPbPb;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict_FlowSkim dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/FlowSkimRun3DataPbPb.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"FlowSkimRun3DataPbPb", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict_FlowSkim",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_FlowSkim_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_FlowSkim_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict_FlowSkim() {
  TriggerDictionaryInitialization_dict_FlowSkim_Impl();
}
