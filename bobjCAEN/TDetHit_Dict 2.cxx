// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TDetHit_Dict
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
#include "TDetHit.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TDetHit(void *p = 0);
   static void *newArray_TDetHit(Long_t size, void *p);
   static void delete_TDetHit(void *p);
   static void deleteArray_TDetHit(void *p);
   static void destruct_TDetHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TDetHit*)
   {
      ::TDetHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TDetHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TDetHit", ::TDetHit::Class_Version(), "TDetHit.hxx", 15,
                  typeid(::TDetHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TDetHit::Dictionary, isa_proxy, 4,
                  sizeof(::TDetHit) );
      instance.SetNew(&new_TDetHit);
      instance.SetNewArray(&newArray_TDetHit);
      instance.SetDelete(&delete_TDetHit);
      instance.SetDeleteArray(&deleteArray_TDetHit);
      instance.SetDestructor(&destruct_TDetHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TDetHit*)
   {
      return GenerateInitInstanceLocal((::TDetHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TDetHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TDetHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TDetHit::Class_Name()
{
   return "TDetHit";
}

//______________________________________________________________________________
const char *TDetHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDetHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TDetHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDetHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TDetHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDetHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TDetHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDetHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TDetHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class TDetHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TDetHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(TDetHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TDetHit(void *p) {
      return  p ? new(p) ::TDetHit : new ::TDetHit;
   }
   static void *newArray_TDetHit(Long_t nElements, void *p) {
      return p ? new(p) ::TDetHit[nElements] : new ::TDetHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_TDetHit(void *p) {
      delete ((::TDetHit*)p);
   }
   static void deleteArray_TDetHit(void *p) {
      delete [] ((::TDetHit*)p);
   }
   static void destruct_TDetHit(void *p) {
      typedef ::TDetHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TDetHit

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 471,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_TDetHit_Dict_Impl() {
    static const char* headers[] = {
"TDetHit.hxx",
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
#line 1 "TDetHit_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TDetHit.hxx")))  TDetHit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TDetHit_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TDetHit.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TDetHit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TDetHit_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TDetHit_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TDetHit_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TDetHit_Dict() {
  TriggerDictionaryInitialization_TDetHit_Dict_Impl();
}
