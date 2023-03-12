'''
https://stackoverflow.com/questions/4801366/convert-rgb-values-to-integer
https://stackoverflow.com/questions/2262100/rgb-int-to-rgb-python
'''
# 色彩编码转换，无需改动

def rgb2int(Red: int, Green: int, Blue: int):
    RGBint = Red
    RGBint = (RGBint << 8) + Green
    RGBint = (RGBint << 8) + Blue
    return RGBint


def int2rgb(RGBint: int):
    Blue = RGBint & 255
    Green = (RGBint >> 8) & 255
    Red = (RGBint >> 16) & 255
    return (Red, Green, Blue)
