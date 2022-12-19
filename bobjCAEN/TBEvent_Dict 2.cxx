// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBEvent_Dict
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
#include "TBEvent.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TBEvent(void *p = 0);
   static void *newArray_TBEvent(Long_t size, void *p);
   static void delete_TBEvent(void *p);
   static void deleteArray_TBEvent(void *p);
   static void destruct_TBEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBEvent*)
   {
      ::TBEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBEvent", ::TBEvent::Class_Version(), "TBEvent.hxx", 18,
                  typeid(::TBEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TBEvent) );
      instance.SetNew(&new_TBEvent);
      instance.SetNewArray(&newArray_TBEvent);
      instance.SetDelete(&delete_TBEvent);
      instance.SetDeleteArray(&deleteArray_TBEvent);
      instance.SetDestructor(&destruct_TBEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBEvent*)
   {
      return GenerateInitInstanceLocal((::TBEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBEvent::Class_Name()
{
   return "TBEvent";
}

//______________________________________________________________________________
const char *TBEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBEvent(void *p) {
      return  p ? new(p) ::TBEvent : new ::TBEvent;
   }
   static void *newArray_TBEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TBEvent[nElements] : new ::TBEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBEvent(void *p) {
      delete ((::TBEvent*)p);
   }
   static void deleteArray_TBEvent(void *p) {
      delete [] ((::TBEvent*)p);
   }
   static void destruct_TBEvent(void *p) {
      typedef ::TBEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBEvent

namespace ROOT {
   static TClass *vectorlETDetHitgR_Dictionary();
   static void vectorlETDetHitgR_TClassManip(TClass*);
   static void *new_vectorlETDetHitgR(void *p = 0);
   static void *newArray_vectorlETDetHitgR(Long_t size, void *p);
   static void delete_vectorlETDetHitgR(void *p);
   static void deleteArray_vectorlETDetHitgR(void *p);
   static void destruct_vectorlETDetHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TDetHit>*)
   {
      vector<TDetHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TDetHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TDetHit>", -2, "vector", 471,
                  typeid(vector<TDetHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETDetHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TDetHit>) );
      instance.SetNew(&new_vectorlETDetHitgR);
      instance.SetNewArray(&newArray_vectorlETDetHitgR);
      instance.SetDelete(&delete_vectorlETDetHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlETDetHitgR);
      instance.SetDestructor(&destruct_vectorlETDetHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TDetHit> >()));

      ::ROOT::AddClassAlternate("vector<TDetHit>","std::__1::vector<TDetHit, std::__1::allocator<TDetHit> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TDetHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETDetHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TDetHit>*)0x0)->GetClass();
      vectorlETDetHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETDetHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETDetHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TDetHit> : new vector<TDetHit>;
   }
   static void *newArray_vectorlETDetHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TDetHit>[nElements] : new vector<TDetHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETDetHitgR(void *p) {
      delete ((vector<TDetHit>*)p);
   }
   static void deleteArray_vectorlETDetHitgR(void *p) {
      delete [] ((vector<TDetHit>*)p);
   }
   static void destruct_vectorlETDetHitgR(void *p) {
      typedef vector<TDetHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TDetHit>

namespace {
  void TriggerDictionaryInitialization_TBEvent_Dict_Impl() {
    static const char* headers[] = {
"TBEvent.hxx",
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
#line 1 "TBEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TDetHit.hxx")))  __attribute__((annotate("$clingAutoload$TBEvent.hxx")))  TDetHit;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TBEvent.hxx")))  TBEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBEvent_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBEvent.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBEvent_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBEvent_Dict() {
  TriggerDictionaryInitialization_TBEvent_Dict_Impl();
}
