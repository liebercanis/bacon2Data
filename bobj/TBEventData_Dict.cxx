// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBEventData_Dict
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

// Header files passed as explicit arguments
#include "TBEventData.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TBEventData(void *p = 0);
   static void *newArray_TBEventData(Long_t size, void *p);
   static void delete_TBEventData(void *p);
   static void deleteArray_TBEventData(void *p);
   static void destruct_TBEventData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBEventData*)
   {
      ::TBEventData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBEventData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBEventData", ::TBEventData::Class_Version(), "TBEventData.hxx", 14,
                  typeid(::TBEventData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBEventData::Dictionary, isa_proxy, 4,
                  sizeof(::TBEventData) );
      instance.SetNew(&new_TBEventData);
      instance.SetNewArray(&newArray_TBEventData);
      instance.SetDelete(&delete_TBEventData);
      instance.SetDeleteArray(&deleteArray_TBEventData);
      instance.SetDestructor(&destruct_TBEventData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBEventData*)
   {
      return GenerateInitInstanceLocal((::TBEventData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBEventData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBEventData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBEventData::Class_Name()
{
   return "TBEventData";
}

//______________________________________________________________________________
const char *TBEventData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBEventData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBEventData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBEventData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBEventData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBEventData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBEventData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBEventData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBEventData::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBEventData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBEventData::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBEventData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBEventData(void *p) {
      return  p ? new(p) ::TBEventData : new ::TBEventData;
   }
   static void *newArray_TBEventData(Long_t nElements, void *p) {
      return p ? new(p) ::TBEventData[nElements] : new ::TBEventData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBEventData(void *p) {
      delete ((::TBEventData*)p);
   }
   static void deleteArray_TBEventData(void *p) {
      delete [] ((::TBEventData*)p);
   }
   static void destruct_TBEventData(void *p) {
      typedef ::TBEventData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBEventData

namespace {
  void TriggerDictionaryInitialization_TBEventData_Dict_Impl() {
    static const char* headers[] = {
"TBEventData.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.24.06/include/",
"/Users/gold/Documents/GitHub/bacon2Data/bobj/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TBEventData_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TBEventData.hxx")))  TBEventData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBEventData_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBEventData.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBEventData", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBEventData_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBEventData_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBEventData_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBEventData_Dict() {
  TriggerDictionaryInitialization_TBEventData_Dict_Impl();
}
