using YAXArrays

using YAXArrays, DimensionalData
axlist = (
    Dim{:time}(range(1, 20, length=20)),
    X(range(1, 10, length=10)),
    Y(range(1, 5, length=15)),
    Dim{:Variable}(["var1", "var2"]))


data = rand(20, 10, 15, 2)

ds = YAXArray(axlist, data, props)