# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright 2016 Christoph Wiedemann

import struct, sys, pprint, glob
import numpy

def read_exr(filename, showinfo = True):
    with open(filename, 'rb') as f:
        magic_number, = struct.unpack("I", f.read(4))
        assert magic_number == 20000630
        ffversion, = struct.unpack("B", f.read(1))
        assert ffversion == 2
        ftype, = struct.unpack("I", b"\0" + f.read(3))
        assert ftype == 0 or (ftype & 2) == 2
        # we assume single part scan line images
        
        # header is a sequence of attributes ended by a null byte
        def read_attr(f):
        
            def read_str(f):
                res = b""
                while 1:
                    c = f.read(1)
                    if c == b"\0":
                        return res.decode('ascii')
                    else:
                        res += c
                        
            def read_channel_list(f):
            
                def read_channel(f):
                    name = read_str(f)
                    if name == "":
                        return None
                    PIXEL_TYPE_ENUM = {0:"UINT", 1:"HALF", 2:"FLOAT"}
                    pixel_type = PIXEL_TYPE_ENUM[struct.unpack("i", f.read(4))[0]]
                    pLinear, = struct.unpack("B", f.read(1))
                    reserved = f.read(3)
                    xSampling, ySampling = struct.unpack("ii", f.read(8))
                    return name, pixel_type, pLinear, xSampling, ySampling
                    
                res = []
                while 1:
                    c = read_channel(f)
                    if c is None:
                        return res
                    res.append(c)
                        
            name = read_str(f)
            if name == "":
                return None
            type = read_str(f)
            size, = struct.unpack("i", f.read(4))
            if type == "box2i":
                value = struct.unpack("iiii", f.read(4*4))
            elif type == "box2f":
                value = struct.unpack("ffff", f.read(4*4))
            elif type == "chlist":
                value = read_channel_list(f)
            elif type == "chromaticities":
                value = struct.unpack("ffffffff", f.read(4*8))
            elif type == "compression":
                COMPRESSION_ENUM = {0:"NONE", 1:"RLE", 2:"ZIPS", 3:"ZIP", 4:"PIZ", 5:"PXR24", 6:"B44", 7:"B44A"}
                value = COMPRESSION_ENUM[struct.unpack("B", f.read(1))[0]]
            elif type == "double":
                value = struct.unpack("d", f.read(8))
            elif type == "envmap":
                ENVMAP_ENUM = {0:"LATLONG", 1:"CUBE"}
                value = ENVMAP_ENUM[struct.unpack("B", f.read(1))[0]]
            elif type == "float":
                value = struct.unpack("f", f.read(4))
            elif type == "int":
                value = struct.unpack("i", f.read(4))
            elif type == "keycode":
                value = struct.unpack("iiiiiii", f.read(4*7))
            elif type == "lineOrder":
                LINE_ORDER_ENUM = {0:"INCREASING_Y", 1:"DECREASING_Y", 2:"RANDOM_Y"}
                value = LINE_ORDER_ENUM[struct.unpack("B", f.read(1))[0]]
            elif type == "m33f":
                value = struct.unpack("fffffffff", f.read(9*4))
            elif type == "m44f":
                value = struct.unpack("ffffffffffffffff", f.read(16*4))
            elif type == "rational":
                value = struct.unpack("iI", f.read(2*4))
            elif type == "string":
                value = read_str(f)
            elif type == "stringvector":
                assert False # not yet implemented
            elif type == "tiledesc":
                value = struct.unpack("IIB", f.read(9))
            elif type == "timecode":
                value = struct.unpack("II", f.read(8))
            elif type == "v2i":
                value = struct.unpack("ii", f.read(8))
            elif type == "v2f":
                value = struct.unpack("ff", f.read(8))
            elif type == "v3i":
                value = struct.unpack("iii", f.read(12))
            elif type == "v3f":
                value = struct.unpack("fff", f.read(12))
            
            return name, type, size, value

        attrs = {}
        while 1:
            a = read_attr(f)
            if a is None:
                break
            attrs[a[0]] = a[-1]
        
        if showinfo:
            pprint.pprint(attrs)
        
        # offset table (scanline assumed)
        comp = attrs["compression"]
        dw = attrs["dataWindow"]
        if comp in ["NONE", "RLE", "ZIPS"]:
            nscan_lines_per_block = 1
        elif comp in ["ZIP", "PXR24"]:
            nscan_lines_per_block = 16
        elif comp in ["PIZ", "B44", "B44A"]:
            nscan_lines_per_block = 32
        cnt_offsets = (dw[3] - dw[1] + 1 + nscan_lines_per_block-1)//nscan_lines_per_block
        
        offsets = struct.unpack("Q"*cnt_offsets, f.read(8*cnt_offsets))
        #print(offsets)
        
        def uncompress_NONE(pixel_data, attr):
            dw = attr["dataWindow"]
            width = dw[2] - dw[0] + 1
            chlist = attr["channels"]
            nch = len(chlist)
            # float32 is assumed
            assert len(pixel_data) == nch * width * 4
            R = numpy.fromstring(pixel_data, dtype='<f4')
            R.shape = [nch, width]
            return R

        try:
            uncompressor = locals()["uncompress_" + comp]
        except KeyError:
            assert False # uncompressor not yet implemented
            
        R = numpy.zeros( (dw[3] - dw[1] + 1, dw[2] - dw[0] + 1, len(attrs['channels']) ), dtype="f")
        for o in offsets:
            f.seek(o)
            ycoord, = struct.unpack("i", f.read(4))
            pxdatasize, = struct.unpack("i", f.read(4))
            pixel_data = f.read(pxdatasize)
            pixels = uncompressor(pixel_data, attrs)
            for y in range(ycoord, ycoord + nscan_lines_per_block):
                for c in range(len(attrs['channels'])):
                    R[y,:,c] = pixels[c,:]
        return R
        
if __name__ == "__main__":
    if "-show" in sys.argv:
        sys.argv.remove("-show")
        show = True
    else:
        show = False

    infns = sys.argv[1]
    infns = glob.glob(infns)
    
    for infn in infns:
        if len(sys.argv) >= 3:
            outfn = sys.argv[2]  
        else:
            outfn = infn
            if outfn[-4:].lower() == ".exr":
                outfn = outfn[:-4]
            outfn += ".mat"
        
        print("Converting", infn, "to", outfn)

        X = read_exr(infn)
        import scipy.io as sio 
        sio.savemat(outfn, {'data':X})
        
        if show:
            import matplotlib.pyplot as plt
            plt.imshow(X[:,:,0] + X[:,:,1] + X[:,:,2])
            plt.colorbar()
            plt.show()
        
        