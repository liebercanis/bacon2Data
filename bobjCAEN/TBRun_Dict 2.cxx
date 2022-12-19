// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBRun_Dict
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
#include "TBRun.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TBRun(void *p = 0);
   static void *newArray_TBRun(Long_t size, void *p);
   static void delete_TBRun(void *p);
   static void deleteArray_TBRun(void *p);
   static void destruct_TBRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBRun*)
   {
      ::TBRun *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBRun >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBRun", ::TBRun::Class_Version(), "TBRun.hxx", 19,
                  typeid(::TBRun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBRun::Dictionary, isa_proxy, 4,
                  sizeof(::TBRun) );
      instance.SetNew(&new_TBRun);
      instance.SetNewArray(&newArray_TBRun);
      instance.SetDelete(&delete_TBRun);
      instance.SetDeleteArray(&deleteArray_TBRun);
      instance.SetDestructor(&destruct_TBRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBRun*)
   {
      return GenerateInitInstanceLocal((::TBRun*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBRun*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBRun::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBRun::Class_Name()
{
   return "TBRun";
}

//______________________________________________________________________________
const char *TBRun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRun*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBRun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRun*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBRun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRun*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBRun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRun*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBRun::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBRun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBRun::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBRun::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBRun(void *p) {
      return  p ? new(p) ::TBRun : new ::TBRun;
   }
   static void *newArray_TBRun(Long_t nElements, void *p) {
      return p ? new(p) ::TBRun[nElements] : new ::TBRun[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBRun(void *p) {
      delete ((::TBRun*)p);
   }
   static void deleteArray_TBRun(void *p) {
      delete [] ((::TBRun*)p);
   }
   static void destruct_TBRun(void *p) {
      typedef ::TBRun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBRun

namespace ROOT {
   static TClass *vectorlETDetmUgR_Dictionary();
   static void vectorlETDetmUgR_TClassManip(TClass*);
   static void *new_vectorlETDetmUgR(void *p = 0);
   static void *newArray_vectorlETDetmUgR(Long_t size, void *p);
   static void delete_vectorlETDetmUgR(void *p);
   static void deleteArray_vectorlETDetmUgR(void *p);
   static void destruct_vectorlETDetmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TDet*>*)
   {
      vector<TDet*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TDet*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TDet*>", -2, "vector", 471,
                  typeid(vector<TDet*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETDetmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TDet*>) );
      instance.SetNew(&new_vectorlETDetmUgR);
      instance.SetNewArray(&newArray_vectorlETDetmUgR);
      instance.SetDelete(&delete_vectorlETDetmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETDetmUgR);
      instance.SetDestructor(&destruct_vectorlETDetmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TDet*> >()));

      ::ROOT::AddClassAlternate("vector<TDet*>","std::__1::vector<TDet*, std::__1::allocator<TDet*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TDet*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETDetmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TDet*>*)0x0)->GetClass();
      vectorlETDetmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETDetmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETDetmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TDet*> : new vector<TDet*>;
   }
   static void *newArray_vectorlETDetmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TDet*>[nElements] : new vector<TDet*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETDetmUgR(void *p) {
      delete ((vector<TDet*>*)p);
   }
   static void deleteArray_vectorlETDetmUgR(void *p) {
      delete [] ((vector<TDet*>*)p);
   }
   static void destruct_vectorlETDetmUgR(void *p) {
      typedef vector<TDet*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TDet*>

namespace {
  void TriggerDictionaryInitialization_TBRun_Dict_Impl() {
    static const char* headers[] = {
"TBRun.hxx",
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
#line 1 "TBRun_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TDet.hxx")))  __attribute__((annotate("$clingAutoload$TBRun.hxx")))  TDet;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TBRun.hxx")))  TBRun;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBRun_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBRun.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBRun", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBRun_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBRun_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBRun_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBRun_Dict() {
  TriggerDictionaryInitialization_TBRun_Dict_Impl();
}
