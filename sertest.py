from fragment_index import serializer
import io

x = serializer.SerializerBlockList()
x.append('foo0' * 12)
x.append('1bar' * 6)
buff = io.BytesIO()
x.write(buff)
print(buff.getvalue())
buff.seek(0)
y = x.read(buff)
print(bytearray(y[0]))

x.to_file("foo.bin")

with open("foo.bin", 'rb') as fh:
    z = x.read(fh)
print(bytearray(y[0]))
