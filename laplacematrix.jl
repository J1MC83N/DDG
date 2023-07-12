# Adapted from Digital Domain's DiscreteDifferentialGeometry.jl


# Copyright 2020 Digital Domain 3.0
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

import Meshes: laplacematrix


function laplacematrix(mesh::IHTriMesh{Dim,T};shift::T=zero(T)) where {Dim,T}
    l = 9*nfaces(mesh) + !iszero(shift)*nvertices(mesh)
    P = mesh.vertices
    indexI = Vector{Int}(undef,l)
    indexJ = Vector{Int}(undef,l)
    element = Vector{T}(undef,l)
    s = 1
    @inbounds for fid in faceids(mesh)
        i,j,k = FVIterator(mesh,fid)
        Pᵢ = P[i]
        Pⱼ = P[j]
        Pₖ = P[k]
        ij = Pⱼ-Pᵢ
        jk = Pₖ-Pⱼ
        ki = Pᵢ-Pₖ
        a = -_cotan(ij,-ki)
        b = -_cotan(jk,-ij)
        c = -_cotan(ki,-jk)
        #L[i,i] += b+c
        indexI[s] = i; indexJ[s] = i; element[s] = b+c; s+=1
        #L[i,j] -= c
        indexI[s] = i; indexJ[s] = j; element[s] = -c; s+=1
        #L[i,k] -= b
        indexI[s] = i; indexJ[s] = k; element[s] = -b; s+=1
        #L[j,i] -= c
        indexI[s] = j; indexJ[s] = i; element[s] = -c; s+=1
        #L[j,j] += c+a
        indexI[s] = j; indexJ[s] = j; element[s] = c+a; s+=1
        #L[j,k] -= a
        indexI[s] = j; indexJ[s] = k; element[s] = -a; s+=1
        #L[k,i] -= b
        indexI[s] = k; indexJ[s] = i; element[s] = -b; s+=1
        #L[k,j] -= a
        indexI[s] = k; indexJ[s] = j; element[s] = -a; s+=1
        #L[k,k] += a+b
        indexI[s] = k; indexJ[s] = k; element[s] = a+b; s+=1
    end
    if !iszero(shift)
        @inbounds for vid in vertexids(mesh)
            indexI[s] = vid; indexJ[s] = vid; element[s] = shift; s+=1
        end
    end
    sparse(indexI, indexJ, element)
end