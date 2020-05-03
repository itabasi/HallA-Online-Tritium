// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME onlineDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "online.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *OnlineConfig_Dictionary();
   static void OnlineConfig_TClassManip(TClass*);
   static void *new_OnlineConfig(void *p = 0);
   static void *newArray_OnlineConfig(Long_t size, void *p);
   static void delete_OnlineConfig(void *p);
   static void deleteArray_OnlineConfig(void *p);
   static void destruct_OnlineConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::OnlineConfig*)
   {
      ::OnlineConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::OnlineConfig));
      static ::ROOT::TGenericClassInfo 
         instance("OnlineConfig", "online.h", 66,
                  typeid(::OnlineConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &OnlineConfig_Dictionary, isa_proxy, 0,
                  sizeof(::OnlineConfig) );
      instance.SetNew(&new_OnlineConfig);
      instance.SetNewArray(&newArray_OnlineConfig);
      instance.SetDelete(&delete_OnlineConfig);
      instance.SetDeleteArray(&deleteArray_OnlineConfig);
      instance.SetDestructor(&destruct_OnlineConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::OnlineConfig*)
   {
      return GenerateInitInstanceLocal((::OnlineConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::OnlineConfig*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *OnlineConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::OnlineConfig*)0x0)->GetClass();
      OnlineConfig_TClassManip(theClass);
   return theClass;
   }

   static void OnlineConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *OnlineGUI_Dictionary();
   static void OnlineGUI_TClassManip(TClass*);
   static void delete_OnlineGUI(void *p);
   static void deleteArray_OnlineGUI(void *p);
   static void destruct_OnlineGUI(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::OnlineGUI*)
   {
      ::OnlineGUI *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::OnlineGUI));
      static ::ROOT::TGenericClassInfo 
         instance("OnlineGUI", "online.h", 113,
                  typeid(::OnlineGUI), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &OnlineGUI_Dictionary, isa_proxy, 0,
                  sizeof(::OnlineGUI) );
      instance.SetDelete(&delete_OnlineGUI);
      instance.SetDeleteArray(&deleteArray_OnlineGUI);
      instance.SetDestructor(&destruct_OnlineGUI);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::OnlineGUI*)
   {
      return GenerateInitInstanceLocal((::OnlineGUI*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::OnlineGUI*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *OnlineGUI_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::OnlineGUI*)0x0)->GetClass();
      OnlineGUI_TClassManip(theClass);
   return theClass;
   }

   static void OnlineGUI_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_OnlineConfig(void *p) {
      return  p ? new(p) ::OnlineConfig : new ::OnlineConfig;
   }
   static void *newArray_OnlineConfig(Long_t nElements, void *p) {
      return p ? new(p) ::OnlineConfig[nElements] : new ::OnlineConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_OnlineConfig(void *p) {
      delete ((::OnlineConfig*)p);
   }
   static void deleteArray_OnlineConfig(void *p) {
      delete [] ((::OnlineConfig*)p);
   }
   static void destruct_OnlineConfig(void *p) {
      typedef ::OnlineConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::OnlineConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_OnlineGUI(void *p) {
      delete ((::OnlineGUI*)p);
   }
   static void deleteArray_OnlineGUI(void *p) {
      delete [] ((::OnlineGUI*)p);
   }
   static void destruct_OnlineGUI(void *p) {
      typedef ::OnlineGUI current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::OnlineGUI

namespace {
  void TriggerDictionaryInitialization_onlineDict_Impl() {
    static const char* headers[] = {
"online.h",
0
    };
    static const char* includePaths[] = {
"/cern/root-6.18.00/include",
"/work/dragon/HallA-Online-Tritium/replay_ole0/onlineGUI64/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "onlineDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$online.h")))  OnlineConfig;
class __attribute__((annotate("$clingAutoload$online.h")))  OnlineGUI;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "onlineDict dictionary payload"

#ifndef DICT
  #define DICT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "online.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"OnlineConfig", payloadCode, "@",
"OnlineGUI", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("onlineDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_onlineDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_onlineDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_onlineDict() {
  TriggerDictionaryInitialization_onlineDict_Impl();
}
