#pragma once

#include <variant>
#include <cstddef>

namespace MultiDimVariant {
   // dummy method returning a variant which allows for one additional type than the given variant
   template <class T, typename...Args>
   constexpr auto extend(const std::variant<Args...> &, const T &) -> std::variant<Args..., T> {
      std::variant<Args..., T> ret;
      return ret;
   }

   // helper structor to hold the declaration of the variant
   template <class TypeHelper, std::size_t N>
   struct MakeMultiDimVariant {
      using variant_type = decltype( extend(MakeMultiDimVariant<TypeHelper, N-1>::m_val, typename TypeHelper::template type<N>()) );
      variant_type m_val;
   };

   // specialisation of above helper structure for N=1 i.e. the variant just allowing T<1>
   template <class TypeHelper>
   struct MakeMultiDimVariant<TypeHelper,1> {
      using type_helper = TypeHelper;
      using variant_type = std::variant< typename TypeHelper::template type<1> >;
      variant_type m_val;
   };
}
