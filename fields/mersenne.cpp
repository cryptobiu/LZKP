//
// Created by moriya on 01/10/17.
// Adjusted by Roee on 04/09/18.
//

#include "mersenne.h"

//#include <cstring>

//
//template <>
//TemplateField<ZpMersenneIntElement>::TemplateField(long fieldParam) {
//
//  this->fieldParam = 2147483647;
//  this->elementSizeInBytes = 4;//round up to the next byte
//  this->elementSizeInBits = 31;
//
//  m_ZERO = new ZpMersenneIntElement(0);
//  m_ONE = new ZpMersenneIntElement(1);
//}
//
//template <>
//TemplateField<ZpMersenneLongElement>::TemplateField(long fieldParam) {
//
//  this->elementSizeInBytes = 8;//round up to the next byte
//  this->elementSizeInBits = 61;
//
//  m_ZERO = new ZpMersenneLongElement(0);
//  m_ONE = new ZpMersenneLongElement(1);
//}




//void copy_byte_array_to_byte_vector(const byte *src, int src_len, std::vector<byte> &target_vector, int beginIndex)
//{
//  if ((int) target_vector.size() < beginIndex + src_len)
//    target_vector.resize(beginIndex + src_len);
//  std::memcpy(target_vector.data() + beginIndex, src, src_len);
//  target_vector.insert(target_vector.begin() + beginIndex, src, &src[src_len]);
//}
//
//template <>
//void TemplateField<ZpMersenneIntElement>::elementToBytes(unsigned char *elemenetInBytes, ZpMersenneIntElement &element) {
//  std::memcpy(elemenetInBytes, (byte*)(&element.elem), 4);
//}
//
//template <>
//void TemplateField<ZpMersenneLongElement>::elementToBytes(unsigned char *elemenetInBytes, ZpMersenneLongElement &element) {
//  std::memcpy(elemenetInBytes, (byte*)(&element.elem), 8);
//}
//
//template <>
//void TemplateField<ZpMersenneIntElement>::elementVectorToByteVector(std::vector<ZpMersenneIntElement> &elementVector, std::vector<byte> &byteVector) {
//  copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size() * elementSizeInBytes, byteVector, 0);
//}
//
//template <>
//void TemplateField<ZpMersenneLongElement>::elementVectorToByteVector(std::vector<ZpMersenneLongElement> &elementVector, std::vector<byte> &byteVector) {
//  copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size() * elementSizeInBytes, byteVector, 0);
//}
//
//template <>
//ZpMersenneIntElement TemplateField<ZpMersenneIntElement>::bytesToElement(unsigned char *elemenetInBytes) {
//  return ZpMersenneIntElement(*(unsigned int *)elemenetInBytes);
//}
//
//template <>
//ZpMersenneLongElement TemplateField<ZpMersenneLongElement>::bytesToElement(unsigned char *elemenetInBytes) {
//  return ZpMersenneLongElement(*(unsigned long *)elemenetInBytes);
//}