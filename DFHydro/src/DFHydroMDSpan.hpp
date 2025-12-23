#ifndef DFMDSPAN_H
#define DFMDSPAN_H

#include<experimental/mdspan>

namespace DFHydro {

// Example: MDSpan<double,  2, 2> view(data())
template<typename T, int...Ds>
using MDSpan = std::experimental::mdspan<T, std::experimental::extents<int,Ds...>>;
}

#endif // DFMDSPAN_H
