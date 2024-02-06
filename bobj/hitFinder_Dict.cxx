// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME hitFinder_Dict
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
#include "hitFinder.hxx"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *hitFinder_Dictionary();
   static void hitFinder_TClassManip(TClass*);
   static void delete_hitFinder(void *p);
   static void deleteArray_hitFinder(void *p);
   static void destruct_hitFinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::hitFinder*)
   {
      ::hitFinder *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::hitFinder));
      static ::ROOT::TGenericClassInfo 
         instance("hitFinder", "hitFinder.hxx", 50,
                  typeid(::hitFinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &hitFinder_Dictionary, isa_proxy, 4,
                  sizeof(::hitFinder) );
      instance.SetDelete(&delete_hitFinder);
      instance.SetDeleteArray(&deleteArray_hitFinder);
      instance.SetDestructor(&destruct_hitFinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::hitFinder*)
   {
      return GenerateInitInstanceLocal((::hitFinder*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::hitFinder*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *hitFinder_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::hitFinder*)nullptr)->GetClass();
      hitFinder_TClassManip(theClass);
   return theClass;
   }

   static void hitFinder_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_hitFinder(void *p) {
      delete ((::hitFinder*)p);
   }
   static void deleteArray_hitFinder(void *p) {
      delete [] ((::hitFinder*)p);
   }
   static void destruct_hitFinder(void *p) {
      typedef ::hitFinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::hitFinder

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "c++/v1/vector", 478,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));

      ::ROOT::AddClassAlternate("vector<unsigned int>","std::__1::vector<unsigned int, std::__1::allocator<unsigned int>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary();
   static void vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);
   static void destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<unsigned int,unsigned int> >*)
   {
      vector<pair<unsigned int,unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<unsigned int,unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<unsigned int,unsigned int> >", -2, "c++/v1/vector", 478,
                  typeid(vector<pair<unsigned int,unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<unsigned int,unsigned int> >) );
      instance.SetNew(&new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<unsigned int,unsigned int> > >()));

      ::ROOT::AddClassAlternate("vector<pair<unsigned int,unsigned int> >","std::__1::vector<std::__1::pair<unsigned int, unsigned int>, std::__1::allocator<std::__1::pair<unsigned int, unsigned int>>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<pair<unsigned int,unsigned int> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<unsigned int,unsigned int> >*)nullptr)->GetClass();
      vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<unsigned int,unsigned int> > : new vector<pair<unsigned int,unsigned int> >;
   }
   static void *newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<unsigned int,unsigned int> >[nElements] : new vector<pair<unsigned int,unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      delete ((vector<pair<unsigned int,unsigned int> >*)p);
   }
   static void deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      delete [] ((vector<pair<unsigned int,unsigned int> >*)p);
   }
   static void destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      typedef vector<pair<unsigned int,unsigned int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<unsigned int,unsigned int> >

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "c++/v1/vector", 478,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::__1::vector<int, std::__1::allocator<int>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)nullptr)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "c++/v1/vector", 478,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
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

namespace ROOT {
   static TClass *vectorlEcomplexlEdoublegRsPgR_Dictionary();
   static void vectorlEcomplexlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEcomplexlEdoublegRsPgR(void *p = nullptr);
   static void *newArray_vectorlEcomplexlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEcomplexlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEcomplexlEdoublegRsPgR(void *p);
   static void destruct_vectorlEcomplexlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<complex<double> >*)
   {
      vector<complex<double> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<complex<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<complex<double> >", -2, "c++/v1/vector", 478,
                  typeid(vector<complex<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEcomplexlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<complex<double> >) );
      instance.SetNew(&new_vectorlEcomplexlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEcomplexlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEcomplexlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEcomplexlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEcomplexlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<complex<double> > >()));

      ::ROOT::AddClassAlternate("vector<complex<double> >","std::__1::vector<std::__1::complex<double>, std::__1::allocator<std::__1::complex<double>>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<complex<double> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEcomplexlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<complex<double> >*)nullptr)->GetClass();
      vectorlEcomplexlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEcomplexlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEcomplexlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<complex<double> > : new vector<complex<double> >;
   }
   static void *newArray_vectorlEcomplexlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<complex<double> >[nElements] : new vector<complex<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEcomplexlEdoublegRsPgR(void *p) {
      delete ((vector<complex<double> >*)p);
   }
   static void deleteArray_vectorlEcomplexlEdoublegRsPgR(void *p) {
      delete [] ((vector<complex<double> >*)p);
   }
   static void destruct_vectorlEcomplexlEdoublegRsPgR(void *p) {
      typedef vector<complex<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<complex<double> >

namespace ROOT {
   static TClass *vectorlETH1DmUgR_Dictionary();
   static void vectorlETH1DmUgR_TClassManip(TClass*);
   static void *new_vectorlETH1DmUgR(void *p = nullptr);
   static void *newArray_vectorlETH1DmUgR(Long_t size, void *p);
   static void delete_vectorlETH1DmUgR(void *p);
   static void deleteArray_vectorlETH1DmUgR(void *p);
   static void destruct_vectorlETH1DmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH1D*>*)
   {
      vector<TH1D*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH1D*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH1D*>", -2, "c++/v1/vector", 478,
                  typeid(vector<TH1D*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH1DmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TH1D*>) );
      instance.SetNew(&new_vectorlETH1DmUgR);
      instance.SetNewArray(&newArray_vectorlETH1DmUgR);
      instance.SetDelete(&delete_vectorlETH1DmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH1DmUgR);
      instance.SetDestructor(&destruct_vectorlETH1DmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH1D*> >()));

      ::ROOT::AddClassAlternate("vector<TH1D*>","std::__1::vector<TH1D*, std::__1::allocator<TH1D*>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TH1D*>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH1DmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TH1D*>*)nullptr)->GetClass();
      vectorlETH1DmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH1DmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH1DmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1D*> : new vector<TH1D*>;
   }
   static void *newArray_vectorlETH1DmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1D*>[nElements] : new vector<TH1D*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH1DmUgR(void *p) {
      delete ((vector<TH1D*>*)p);
   }
   static void deleteArray_vectorlETH1DmUgR(void *p) {
      delete [] ((vector<TH1D*>*)p);
   }
   static void destruct_vectorlETH1DmUgR(void *p) {
      typedef vector<TH1D*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TH1D*>

namespace ROOT {
   static TClass *maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = nullptr);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "c++/v1/map", 912,
                  typeid(map<int,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));

      ::ROOT::AddClassAlternate("map<int,int>","std::__1::map<int, int, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, int>>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)nullptr)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

namespace ROOT {
   static TClass *maplEdoublecOTDetHitgR_Dictionary();
   static void maplEdoublecOTDetHitgR_TClassManip(TClass*);
   static void *new_maplEdoublecOTDetHitgR(void *p = nullptr);
   static void *newArray_maplEdoublecOTDetHitgR(Long_t size, void *p);
   static void delete_maplEdoublecOTDetHitgR(void *p);
   static void deleteArray_maplEdoublecOTDetHitgR(void *p);
   static void destruct_maplEdoublecOTDetHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<double,TDetHit>*)
   {
      map<double,TDetHit> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<double,TDetHit>));
      static ::ROOT::TGenericClassInfo 
         instance("map<double,TDetHit>", -2, "c++/v1/map", 912,
                  typeid(map<double,TDetHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEdoublecOTDetHitgR_Dictionary, isa_proxy, 0,
                  sizeof(map<double,TDetHit>) );
      instance.SetNew(&new_maplEdoublecOTDetHitgR);
      instance.SetNewArray(&newArray_maplEdoublecOTDetHitgR);
      instance.SetDelete(&delete_maplEdoublecOTDetHitgR);
      instance.SetDeleteArray(&deleteArray_maplEdoublecOTDetHitgR);
      instance.SetDestructor(&destruct_maplEdoublecOTDetHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<double,TDetHit> >()));

      ::ROOT::AddClassAlternate("map<double,TDetHit>","std::__1::map<double, TDetHit, std::__1::less<double>, std::__1::allocator<std::__1::pair<double const, TDetHit>>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<double,TDetHit>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEdoublecOTDetHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<double,TDetHit>*)nullptr)->GetClass();
      maplEdoublecOTDetHitgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEdoublecOTDetHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEdoublecOTDetHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<double,TDetHit> : new map<double,TDetHit>;
   }
   static void *newArray_maplEdoublecOTDetHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<double,TDetHit>[nElements] : new map<double,TDetHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEdoublecOTDetHitgR(void *p) {
      delete ((map<double,TDetHit>*)p);
   }
   static void deleteArray_maplEdoublecOTDetHitgR(void *p) {
      delete [] ((map<double,TDetHit>*)p);
   }
   static void destruct_maplEdoublecOTDetHitgR(void *p) {
      typedef map<double,TDetHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<double,TDetHit>

namespace {
  void TriggerDictionaryInitialization_hitFinder_Dict_Impl() {
    static const char* headers[] = {
"hitFinder.hxx",
nullptr
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/usr/local/root-v6-26-10-buildFFTNew/include/",
"/Users/gold/GOOGLE2/bacon2Data/bobj/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "hitFinder_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$hitFinder.hxx")))  hitFinder;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "hitFinder_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "hitFinder.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"hitFinder", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("hitFinder_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_hitFinder_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_hitFinder_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_hitFinder_Dict() {
  TriggerDictionaryInitialization_hitFinder_Dict_Impl();
}
