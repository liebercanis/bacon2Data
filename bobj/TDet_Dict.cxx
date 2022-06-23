// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TDet_Dict
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
#include "TDet.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TDet(void *p = 0);
   static void *newArray_TDet(Long_t size, void *p);
   static void delete_TDet(void *p);
   static void deleteArray_TDet(void *p);
   static void destruct_TDet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TDet*)
   {
      ::TDet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TDet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TDet", ::TDet::Class_Version(), "TDet.hxx", 16,
                  typeid(::TDet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TDet::Dictionary, isa_proxy, 4,
                  sizeof(::TDet) );
      instance.SetNew(&new_TDet);
      instance.SetNewArray(&newArray_TDet);
      instance.SetDelete(&delete_TDet);
      instance.SetDeleteArray(&deleteArray_TDet);
      instance.SetDestructor(&destruct_TDet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TDet*)
   {
      return GenerateInitInstanceLocal((::TDet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TDet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TDet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TDet::Class_Name()
{
   return "TDet";
}

//______________________________________________________________________________
const char *TDet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TDet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TDet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TDet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TDet::Streamer(TBuffer &R__b)
{
   // Stream an object of class TDet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TDet::Class(),this);
   } else {
      R__b.WriteClassBuffer(TDet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TDet(void *p) {
      return  p ? new(p) ::TDet : new ::TDet;
   }
   static void *newArray_TDet(Long_t nElements, void *p) {
      return p ? new(p) ::TDet[nElements] : new ::TDet[nElements];
   }
   // Wrapper around operator delete
   static void delete_TDet(void *p) {
      delete ((::TDet*)p);
   }
   static void deleteArray_TDet(void *p) {
      delete [] ((::TDet*)p);
   }
   static void destruct_TDet(void *p) {
      typedef ::TDet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TDet

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
  void TriggerDictionaryInitialization_TDet_Dict_Impl() {
    static const char* headers[] = {
"TDet.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.24.06/include/",
"/Users/gold/Documents/GitHub/bacon2Data/bacon2Data/bobj/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TDet_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TDetHit.hxx")))  __attribute__((annotate("$clingAutoload$TDet.hxx")))  TDetHit;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TDet.hxx")))  TDet;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TDet_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TDet.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TDet", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TDet_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TDet_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TDet_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TDet_Dict() {
  TriggerDictionaryInitialization_TDet_Dict_Impl();
}
