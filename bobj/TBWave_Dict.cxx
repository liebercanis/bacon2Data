// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBWave_Dict
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
#include "TBWave.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TBWave(void *p = 0);
   static void *newArray_TBWave(Long_t size, void *p);
   static void delete_TBWave(void *p);
   static void deleteArray_TBWave(void *p);
   static void destruct_TBWave(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBWave*)
   {
      ::TBWave *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBWave >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBWave", ::TBWave::Class_Version(), "TBWave.hxx", 17,
                  typeid(::TBWave), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBWave::Dictionary, isa_proxy, 4,
                  sizeof(::TBWave) );
      instance.SetNew(&new_TBWave);
      instance.SetNewArray(&newArray_TBWave);
      instance.SetDelete(&delete_TBWave);
      instance.SetDeleteArray(&deleteArray_TBWave);
      instance.SetDestructor(&destruct_TBWave);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBWave*)
   {
      return GenerateInitInstanceLocal((::TBWave*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBWave*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBWave::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBWave::Class_Name()
{
   return "TBWave";
}

//______________________________________________________________________________
const char *TBWave::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBWave*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBWave::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBWave*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBWave::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBWave*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBWave::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBWave*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBWave::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBWave.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBWave::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBWave::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBWave(void *p) {
      return  p ? new(p) ::TBWave : new ::TBWave;
   }
   static void *newArray_TBWave(Long_t nElements, void *p) {
      return p ? new(p) ::TBWave[nElements] : new ::TBWave[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBWave(void *p) {
      delete ((::TBWave*)p);
   }
   static void deleteArray_TBWave(void *p) {
      delete [] ((::TBWave*)p);
   }
   static void destruct_TBWave(void *p) {
      typedef ::TBWave current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBWave

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
  void TriggerDictionaryInitialization_TBWave_Dict_Impl() {
    static const char* headers[] = {
"TBWave.hxx",
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
#line 1 "TBWave_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TBWave.hxx")))  TBWave;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBWave_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBWave.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBWave", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBWave_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBWave_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBWave_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBWave_Dict() {
  TriggerDictionaryInitialization_TBWave_Dict_Impl();
}
