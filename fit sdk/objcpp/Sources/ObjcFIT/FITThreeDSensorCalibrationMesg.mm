////////////////////////////////////////////////////////////////////////////////
// The following FIT Protocol software provided may be used with FIT protocol
// devices only and remains the copyrighted property of Garmin Canada Inc.
// The software is being provided on an "as-is" basis and as an accommodation,
// and therefore all warranties, representations, or guarantees of any kind
// (whether express, implied or statutory) including, without limitation,
// warranties of merchantability, non-infringement, or fitness for a particular
// purpose, are specifically disclaimed.
//
// Copyright 2021 Garmin International, Inc.
////////////////////////////////////////////////////////////////////////////////
// ****WARNING****  This file is auto-generated!  Do NOT edit this file.
// Profile Version = 21.67Release
// Tag = production/akw/21.67.00-0-gd790f76b
////////////////////////////////////////////////////////////////////////////////


#import "FITMessage+Internal.h"


#import "FITThreeDSensorCalibrationMesg.h"

@implementation FITThreeDSensorCalibrationMesg

- (instancetype)init {
    self = [super initWithFitMesgIndex:fit::Profile::MESG_THREE_D_SENSOR_CALIBRATION];

    return self;
}

// Timestamp 
- (BOOL)isTimestampValid {
	const fit::Field* field = [super getField:253];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid() == FIT_TRUE ? TRUE : FALSE;
}

- (FITDate *)getTimestamp {
    return FITDateFromTimestamp([super getFieldUINT32ValueForField:253 forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setTimestamp:(FITDate *)timestamp {
    [super setFieldUINT32ValueForField:253 andValue:TimestampFromFITDate(timestamp) forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

// SensorType 
- (BOOL)isSensorTypeValid {
	const fit::Field* field = [super getField:0];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid() == FIT_TRUE ? TRUE : FALSE;
}

- (FITSensorType)getSensorType {
    return ([super getFieldENUMValueForField:0 forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setSensorType:(FITSensorType)sensorType {
    [super setFieldENUMValueForField:0 andValue:(sensorType) forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

// CalibrationFactor 
- (BOOL)isCalibrationFactorValid {
	const fit::Field* field = [super getField:1];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid() == FIT_TRUE ? TRUE : FALSE;
}

- (FITUInt32)getCalibrationFactor {
    return ([super getFieldUINT32ValueForField:1 forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setCalibrationFactor:(FITUInt32)calibrationFactor {
    [super setFieldUINT32ValueForField:1 andValue:(calibrationFactor) forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 
// CalibrationFactor - Sub Fields
// AccelCalFactor - CalibrationFactor Field - Sub Field 
- (BOOL)isAccelCalFactorValid
{
    const fit::Field* field = [super getField:1];
    if( FIT_NULL == field ) {
        return FIT_FALSE;
    }

    if(![super canField:1 supportSubField:(FITUInt16)FITThreeDSensorCalibrationMesgCalibrationFactorFieldAccelCalFactorSubField]) {
        return FIT_FALSE;
    }

    return field->IsValueValid(0, FITThreeDSensorCalibrationMesgCalibrationFactorFieldAccelCalFactorSubField) == FIT_TRUE ? TRUE : FALSE;
}

- (FITUInt32)getAccelCalFactor
{
    return ([super getFieldUINT32ValueForField:1 forIndex:0 andSubFieldIndex:FITThreeDSensorCalibrationMesgCalibrationFactorFieldAccelCalFactorSubField]);
}

- (void)setAccelCalFactor:(FITUInt32)accelCalFactor
{
    [super setFieldUINT32ValueForField:1 andValue:(accelCalFactor) forIndex:0 andSubFieldIndex:FITThreeDSensorCalibrationMesgCalibrationFactorFieldAccelCalFactorSubField];
} 
// GyroCalFactor - CalibrationFactor Field - Sub Field 
- (BOOL)isGyroCalFactorValid
{
    const fit::Field* field = [super getField:1];
    if( FIT_NULL == field ) {
        return FIT_FALSE;
    }

    if(![super canField:1 supportSubField:(FITUInt16)FITThreeDSensorCalibrationMesgCalibrationFactorFieldGyroCalFactorSubField]) {
        return FIT_FALSE;
    }

    return field->IsValueValid(0, FITThreeDSensorCalibrationMesgCalibrationFactorFieldGyroCalFactorSubField) == FIT_TRUE ? TRUE : FALSE;
}

- (FITUInt32)getGyroCalFactor
{
    return ([super getFieldUINT32ValueForField:1 forIndex:0 andSubFieldIndex:FITThreeDSensorCalibrationMesgCalibrationFactorFieldGyroCalFactorSubField]);
}

- (void)setGyroCalFactor:(FITUInt32)gyroCalFactor
{
    [super setFieldUINT32ValueForField:1 andValue:(gyroCalFactor) forIndex:0 andSubFieldIndex:FITThreeDSensorCalibrationMesgCalibrationFactorFieldGyroCalFactorSubField];
} 

// CalibrationDivisor 
- (BOOL)isCalibrationDivisorValid {
	const fit::Field* field = [super getField:2];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid() == FIT_TRUE ? TRUE : FALSE;
}

- (FITUInt32)getCalibrationDivisor {
    return ([super getFieldUINT32ValueForField:2 forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setCalibrationDivisor:(FITUInt32)calibrationDivisor {
    [super setFieldUINT32ValueForField:2 andValue:(calibrationDivisor) forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

// LevelShift 
- (BOOL)isLevelShiftValid {
	const fit::Field* field = [super getField:3];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid() == FIT_TRUE ? TRUE : FALSE;
}

- (FITUInt32)getLevelShift {
    return ([super getFieldUINT32ValueForField:3 forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setLevelShift:(FITUInt32)levelShift {
    [super setFieldUINT32ValueForField:3 andValue:(levelShift) forIndex:0 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

// OffsetCal 
- (uint8_t)numOffsetCalValues {
    return [super getFieldNumValuesForField:4 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
}

- (BOOL)isOffsetCalValidforIndex:(uint8_t)index {
	const fit::Field* field = [super getField:4];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid(index) == FIT_TRUE ? TRUE : FALSE;
}

- (FITSInt32)getOffsetCalforIndex:(uint8_t)index {
    return ([super getFieldSINT32ValueForField:4 forIndex:index andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setOffsetCal:(FITSInt32)offsetCal forIndex:(uint8_t)index {
    [super setFieldSINT32ValueForField:4 andValue:(offsetCal) forIndex:index andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

// OrientationMatrix 
- (uint8_t)numOrientationMatrixValues {
    return [super getFieldNumValuesForField:5 andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
}

- (BOOL)isOrientationMatrixValidforIndex:(uint8_t)index {
	const fit::Field* field = [super getField:5];
	if( FIT_NULL == field ) {
		return FALSE;
	}

	return field->IsValueValid(index) == FIT_TRUE ? TRUE : FALSE;
}

- (FITFloat32)getOrientationMatrixforIndex:(uint8_t)index {
    return ([super getFieldFLOAT32ValueForField:5 forIndex:index andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD]);
}

- (void)setOrientationMatrix:(FITFloat32)orientationMatrix forIndex:(uint8_t)index {
    [super setFieldFLOAT32ValueForField:5 andValue:(orientationMatrix) forIndex:index andSubFieldIndex:FIT_SUBFIELD_INDEX_MAIN_FIELD];
} 

@end
