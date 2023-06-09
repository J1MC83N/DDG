# Copyright 2019 Digital Domain 3.0
#
# Licensed under the Apache License, Version 2.0 (the "Apache License")
# with the following modification; you may not use this file except in
# compliance with the Apache License and the following modification to it:
# Section 6. Trademarks. is deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor
#    and its affiliates, except as required to comply with Section 4(c) of
#    the License and to reproduce the content of the NOTICE file.
#
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Apache License with the above modification is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the Apache License for the specific
# language governing permissions and limitations under the Apache License.

"""
A lightweight (primitive), integer-like handle to a resource.  Useful for dispatch.
i.e. pretty much an integer, but with a unique dispatch signature.
"""

export 
Handle, 
@HandleType,
@HandleType_alias

abstract type Handle <: Unsigned end

"""
    @HandleType(n)

Create a new type named n to be used as a handle to your resource.

Example:
julia> @HandleType(Foo)
Foo

julia> foo(x::Integer) = x*x
foo (generic function with 1 method)

julia> foo(Foo(2))
4

julia> foo(x::Foo) = x-x
foo (generic function with 2 methods)

julia> foo(Foo(2))
0
"""
macro HandleType(tname)
  esc(
      quote 
        primitive type $tname <: Handle sizeof(Int)*8 end
        $tname(x::T) where {T<:Integer} = reinterpret($tname, UInt64(unsigned(x)))
      end
     )
end

const _Handle2Alias = Dict{Symbol, Symbol}()
macro HandleType_alias(tname,alias)
  _Handle2Alias[tname] = alias
  esc(
      quote
        primitive type $tname <: Handle sizeof(Int)*8 end
        $tname(x::T) where {T<:Integer} = reinterpret($tname, UInt64(unsigned(x)))
        const $alias = $tname
      end
     )
end

Base.show(io::IO, x::Handle) = print(io, Int(x))
Base.string(x::Handle) = string(Int(x))
Base.ndigits0zpb( c::C, i::Integer ) where C<:Handle = ndigits0zpb( UInt64(c), i )
Base.promote_rule( ::Type{C}, ::Type{I} ) where {C<:Handle,I<:Integer} = C
Base.convert( ::Type{C}, number::Number ) where {C<:Handle} = C(UInt64(number))
Base.typemax(::Type{C}) where {C<:Handle} = C(typemax(UInt64))
Base.typemin(::Type{C}) where {C<:Handle} = C(typemin(UInt64))

Base.Int(x::Handle) = reinterpret(Int,x)
Base.UInt(x::Handle) = reinterpret(UInt,x)
Base.Int32(x::Handle) = Int32(reinterpret(Int,x))
Base.UInt32(x::Handle) = UInt32(reinterpret(UInt,x))
Base.Int128(x::Handle) = Int128(reinterpret(Int,x))
Base.UInt128(x::Handle) = UInt128(reinterpret(UInt,x))

Base.hash(d::Handle, x::UInt64) = hash(reinterpret(UInt64,d),x)

for op in (:+, :-, :*, :/, :mod, :div, :rem, :max, :min, :xor)
  @eval Base.$op(a::H, b::H) where {H<:Handle} = H($op(UInt64(a),UInt64(b)))
end

Base.:(==)(a::Handle, b::Handle) = a === b
for op in (:<, :>, :<=, :>=)
  @eval Base.$op(a::H, b::H) where {H<:Handle} = $op(UInt64(a),UInt64(b))
end

for op in (:<<, :>>), IS in (Int64, UInt64)
  @eval Base.$op(a::H, b::I) where {H<:Handle, I<:$IS} = H($op(UInt64(a),UInt64(b)))
end

for from in Base.BitInteger_types
  if sizeof(UInt64) < sizeof(from)
    @eval Base.rem(x::($from), ::Type{T}) where T<:Handle = Base.trunc_int(T, x)
  end
end
Base.trunc(::Type{H}, x::Real) where H<:Handle = H(trunc(UInt64, x))

Base.parse(::Type{H},str::AbstractString) where H<:Handle = H(parse(Int,str))

Base.sub_with_overflow( a::H, b::H ) where {H<:Handle} = 
  ((r,f)->(H(r),f))(Base.sub_with_overflow(UInt64(a), UInt64(b))...)
Base.add_with_overflow( a::H, b::H ) where {H<:Handle} = 
  ((r,f)->(H(r),f))(Base.add_with_overflow(UInt64(a), UInt64(b))...)

import Random
Random.uint_sup(::Type{H}) where H<:Handle = H
Random.Sampler(::Type{<:Random.AbstractRNG}, r::AbstractUnitRange{H},
        ::Random.Repetition) where H<:Handle = Random.SamplerRangeNDL(r)
Base.widen(::Type{H}) where H<:Handle = UInt128
Base.rand(rng::Random.AbstractRNG, ::Random.SamplerType{H}) where H<:Handle = H(rand(rng,Random.SamplerType{UInt64}()))
Base.promote_rule(::Type{H},::Type{Int128}) where H<:Handle = Int128
Base.promote_rule(::Type{H},::Type{UInt128}) where H<:Handle = UInt128