// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TBRawEvent_Dict
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
#include "TBRawEvent.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TBRawEvent(void *p = 0);
   static void *newArray_TBRawEvent(Long_t size, void *p);
   static void delete_TBRawEvent(void *p);
   static void deleteArray_TBRawEvent(void *p);
   static void destruct_TBRawEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBRawEvent*)
   {
      ::TBRawEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBRawEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBRawEvent", ::TBRawEvent::Class_Version(), "TBRawEvent.hxx", 15,
                  typeid(::TBRawEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBRawEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TBRawEvent) );
      instance.SetNew(&new_TBRawEvent);
      instance.SetNewArray(&newArray_TBRawEvent);
      instance.SetDelete(&delete_TBRawEvent);
      instance.SetDeleteArray(&deleteArray_TBRawEvent);
      instance.SetDestructor(&destruct_TBRawEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBRawEvent*)
   {
      return GenerateInitInstanceLocal((::TBRawEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TBRawEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TBRawEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBRawEvent::Class_Name()
{
   return "TBRawEvent";
}

//______________________________________________________________________________
const char *TBRawEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRawEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBRawEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBRawEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBRawEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRawEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBRawEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBRawEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TBRawEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBRawEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBRawEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBRawEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBRawEvent(void *p) {
      return  p ? new(p) ::TBRawEvent : new ::TBRawEvent;
   }
   static void *newArray_TBRawEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TBRawEvent[nElements] : new ::TBRawEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBRawEvent(void *p) {
      delete ((::TBRawEvent*)p);
   }
   static void deleteArray_TBRawEvent(void *p) {
      delete [] ((::TBRawEvent*)p);
   }
   static void destruct_TBRawEvent(void *p) {
      typedef ::TBRawEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBRawEvent

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 471,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));

      ::ROOT::AddClassAlternate("vector<unsigned short>","std::__1::vector<unsigned short, std::__1::allocator<unsigned short> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned short>*)0x0)->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete ((vector<unsigned short>*)p);
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] ((vector<unsigned short>*)p);
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

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
  void TriggerDictionaryInitialization_TBRawEvent_Dict_Impl() {
    static const char* headers[] = {
"TBRawEvent.hxx",
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
#line 1 "TBRawEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TBRawEvent.hxx")))  TBRawEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TBRawEvent_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TBRawEvent.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TBRawEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TBRawEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TBRawEvent_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TBRawEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TBRawEvent_Dict() {
  TriggerDictionaryInitialization_TBRawEvent_Dict_Impl();
}
