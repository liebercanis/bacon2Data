// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBRawRun_Dict
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
#include "TBRawRun.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void delete_TBRawRun(void *p);
   static void deleteArray_TBRawRun(void *p);
   static void destruct_TBRawRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBRawRun*)
   {
      ::TBRawRun *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBRawRun >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBRawRun", ::TBRawRun::Class_Version(), "TBRawRun.hxx", 17,
                  typeid(::TBRawRun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBRawRun::Dictionary, isa_proxy, 4,
                  sizeof(::TBRawRun) );
      instance.SetDelete(&delete_TBRawRun);
      instance.SetDeleteArray(&deleteArray_TBRawRun);
      instance.SetDestructor(&destruct_TBRawRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBRawRun*)
   {
      return GenerateInitInstanceLocal((::TBRawRun*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBRawRun*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBRawRun::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBRawRun::Class_Name()
{
   return "TBRawRun";
}

//______________________________________________________________________________
const char *TBRawRun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRawRun*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBRawRun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRawRun*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBRawRun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRawRun*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBRawRun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRawRun*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBRawRun::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBRawRun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBRawRun::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBRawRun::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TBRawRun(void *p) {
      delete ((::TBRawRun*)p);
   }
   static void deleteArray_TBRawRun(void *p) {
      delete [] ((::TBRawRun*)p);
   }
   static void destruct_TBRawRun(void *p) {
      typedef ::TBRawRun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBRawRun

namespace ROOT {
   static TClass *vectorlETBRawEventmUgR_Dictionary();
   static void vectorlETBRawEventmUgR_TClassManip(TClass*);
   static void *new_vectorlETBRawEventmUgR(void *p = 0);
   static void *newArray_vectorlETBRawEventmUgR(Long_t size, void *p);
   static void delete_vectorlETBRawEventmUgR(void *p);
   static void deleteArray_vectorlETBRawEventmUgR(void *p);
   static void destruct_vectorlETBRawEventmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TBRawEvent*>*)
   {
      vector<TBRawEvent*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TBRawEvent*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TBRawEvent*>", -2, "vector", 471,
                  typeid(vector<TBRawEvent*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETBRawEventmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TBRawEvent*>) );
      instance.SetNew(&new_vectorlETBRawEventmUgR);
      instance.SetNewArray(&newArray_vectorlETBRawEventmUgR);
      instance.SetDelete(&delete_vectorlETBRawEventmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETBRawEventmUgR);
      instance.SetDestructor(&destruct_vectorlETBRawEventmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TBRawEvent*> >()));

      ::ROOT::AddClassAlternate("vector<TBRawEvent*>","std::__1::vector<TBRawEvent*, std::__1::allocator<TBRawEvent*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TBRawEvent*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETBRawEventmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TBRawEvent*>*)0x0)->GetClass();
      vectorlETBRawEventmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETBRawEventmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETBRawEventmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TBRawEvent*> : new vector<TBRawEvent*>;
   }
   static void *newArray_vectorlETBRawEventmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TBRawEvent*>[nElements] : new vector<TBRawEvent*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETBRawEventmUgR(void *p) {
      delete ((vector<TBRawEvent*>*)p);
   }
   static void deleteArray_vectorlETBRawEventmUgR(void *p) {
      delete [] ((vector<TBRawEvent*>*)p);
   }
   static void destruct_vectorlETBRawEventmUgR(void *p) {
      typedef vector<TBRawEvent*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TBRawEvent*>

namespace {
  void TriggerDictionaryInitialization_TBRawRun_Dict_Impl() {
    static const char* headers[] = {
"TBRawRun.hxx",
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
#line 1 "TBRawRun_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TBRawEvent.hxx")))  __attribute__((annotate("$clingAutoload$TBRawRun.hxx")))  TBRawEvent;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TBRawRun.hxx")))  TBRawRun;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBRawRun_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBRawRun.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBRawRun", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBRawRun_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBRawRun_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBRawRun_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBRawRun_Dict() {
  TriggerDictionaryInitialization_TBRawRun_Dict_Impl();
}
