/**
 * A few useful templates, e.g. transform_if (combination of std::copy_if and std::transform)
 * and LessOf/MoreOf, which turn a unary function into a comparator binary function.
 */
namespace TTWAnalysis {

  template<typename _Arg,typename _UnaryFunction>
  struct _LessOf {
    public:
      _LessOf(_UnaryFunction fun) : m_fun(fun) {}
      bool operator() ( const _Arg& __a, const _Arg& __b ) const
      {
        return m_fun(__a) < m_fun(__b);
      }
    private:
      _UnaryFunction m_fun;
  };
  // factory method
  /* LessOf<_Arg>(_UnaryFunction fun) will comare arguments (_Arg a, _Arg b) by first applying
   * fun, i.e. it will return 'fun(a) < fun(b)'
   */
  template<typename _Arg,typename _UnaryFunction>
  _LessOf<_Arg,_UnaryFunction> LessOf( _UnaryFunction __fun )
  { return _LessOf<_Arg,_UnaryFunction>(std::forward<_UnaryFunction>(__fun)); }

  template<typename _Arg,typename _UnaryFunction>
  struct _MoreOf {
    public:
      _MoreOf(_UnaryFunction fun) : m_fun(fun) {}
      bool operator() ( const _Arg& __a, const _Arg& __b ) const
      {
        return m_fun(__a) > m_fun(__b);
      }
    private:
      _UnaryFunction m_fun;
  };
  // factory method
  /* MoreOf<_Arg>(_UnaryFunction fun) will comare arguments (_Arg a, _Arg b) by first applying
   * fun, i.e. it will return 'fun(a) > fun(b)'
   */
  template<typename _Arg,typename _UnaryFunction>
  _MoreOf<_Arg,_UnaryFunction> MoreOf( _UnaryFunction __fun )
  { return _MoreOf<_Arg,_UnaryFunction>(std::forward<_UnaryFunction>(__fun)); }


  /* workaround the non-composibility of std::copy_if and std::transform
   * (can be removed when using ranges)
   *
   * based on the implementations in libstdc++
   */
  template<typename _InputIterator, typename _OutputIterator,
           typename _Predicate, typename _UnaryOperator>
  _OutputIterator transform_if(_InputIterator __first, _InputIterator __last,
                               _OutputIterator __result,
                               _Predicate __pred,
                               _UnaryOperator __unary_op)
  {
    for (; __first != __last; ++__first)
      if (__pred(*__first))
        {
          *__result = __unary_op(*__first);
          ++__result;
        }
    return __result;
  }
}
