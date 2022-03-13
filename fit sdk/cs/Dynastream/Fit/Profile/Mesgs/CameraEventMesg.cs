#region Copyright
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

#endregion

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using System.IO;
using System.Linq;

namespace Dynastream.Fit
{
    /// <summary>
    /// Implements the CameraEvent profile message.
    /// </summary>
    public class CameraEventMesg : Mesg
    {
        #region Fields
        #endregion

        /// <summary>
        /// Field Numbers for <see cref="CameraEventMesg"/>
        /// </summary>
        public sealed class FieldDefNum
        {
            public const byte Timestamp = 253;
            public const byte TimestampMs = 0;
            public const byte CameraEventType = 1;
            public const byte CameraFileUuid = 2;
            public const byte CameraOrientation = 3;
            public const byte Invalid = Fit.FieldNumInvalid;
        }

        #region Constructors
        public CameraEventMesg() : base(Profile.GetMesg(MesgNum.CameraEvent))
        {
        }

        public CameraEventMesg(Mesg mesg) : base(mesg)
        {
        }
        #endregion // Constructors

        #region Methods
        ///<summary>
        /// Retrieves the Timestamp field
        /// Units: s
        /// Comment: Whole second part of the timestamp.</summary>
        /// <returns>Returns DateTime representing the Timestamp field</returns>
        public DateTime GetTimestamp()
        {
            Object val = GetFieldValue(253, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return TimestampToDateTime(Convert.ToUInt32(val));
            
        }

        /// <summary>
        /// Set Timestamp field
        /// Units: s
        /// Comment: Whole second part of the timestamp.</summary>
        /// <param name="timestamp_">Nullable field value to be set</param>
        public void SetTimestamp(DateTime timestamp_)
        {
            SetFieldValue(253, 0, timestamp_.GetTimeStamp(), Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the TimestampMs field
        /// Units: ms
        /// Comment: Millisecond part of the timestamp.</summary>
        /// <returns>Returns nullable ushort representing the TimestampMs field</returns>
        public ushort? GetTimestampMs()
        {
            Object val = GetFieldValue(0, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToUInt16(val));
            
        }

        /// <summary>
        /// Set TimestampMs field
        /// Units: ms
        /// Comment: Millisecond part of the timestamp.</summary>
        /// <param name="timestampMs_">Nullable field value to be set</param>
        public void SetTimestampMs(ushort? timestampMs_)
        {
            SetFieldValue(0, 0, timestampMs_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the CameraEventType field</summary>
        /// <returns>Returns nullable CameraEventType enum representing the CameraEventType field</returns>
        public CameraEventType? GetCameraEventType()
        {
            object obj = GetFieldValue(1, 0, Fit.SubfieldIndexMainField);
            CameraEventType? value = obj == null ? (CameraEventType?)null : (CameraEventType)obj;
            return value;
        }

        /// <summary>
        /// Set CameraEventType field</summary>
        /// <param name="cameraEventType_">Nullable field value to be set</param>
        public void SetCameraEventType(CameraEventType? cameraEventType_)
        {
            SetFieldValue(1, 0, cameraEventType_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the CameraFileUuid field</summary>
        /// <returns>Returns byte[] representing the CameraFileUuid field</returns>
        public byte[] GetCameraFileUuid()
        {
            byte[] data = (byte[])GetFieldValue(2, 0, Fit.SubfieldIndexMainField);
            return data.Take(data.Length - 1).ToArray();
        }

        ///<summary>
        /// Retrieves the CameraFileUuid field</summary>
        /// <returns>Returns String representing the CameraFileUuid field</returns>
        public String GetCameraFileUuidAsString()
        {
            byte[] data = (byte[])GetFieldValue(2, 0, Fit.SubfieldIndexMainField);
            return data != null ? Encoding.UTF8.GetString(data, 0, data.Length - 1) : null;
        }

        ///<summary>
        /// Set CameraFileUuid field</summary>
        /// <param name="cameraFileUuid_"> field value to be set</param>
        public void SetCameraFileUuid(String cameraFileUuid_)
        {
            byte[] data = Encoding.UTF8.GetBytes(cameraFileUuid_);
            byte[] zdata = new byte[data.Length + 1];
            data.CopyTo(zdata, 0);
            SetFieldValue(2, 0, zdata, Fit.SubfieldIndexMainField);
        }

        
        /// <summary>
        /// Set CameraFileUuid field</summary>
        /// <param name="cameraFileUuid_">field value to be set</param>
        public void SetCameraFileUuid(byte[] cameraFileUuid_)
        {
            SetFieldValue(2, 0, cameraFileUuid_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the CameraOrientation field</summary>
        /// <returns>Returns nullable CameraOrientationType enum representing the CameraOrientation field</returns>
        public CameraOrientationType? GetCameraOrientation()
        {
            object obj = GetFieldValue(3, 0, Fit.SubfieldIndexMainField);
            CameraOrientationType? value = obj == null ? (CameraOrientationType?)null : (CameraOrientationType)obj;
            return value;
        }

        /// <summary>
        /// Set CameraOrientation field</summary>
        /// <param name="cameraOrientation_">Nullable field value to be set</param>
        public void SetCameraOrientation(CameraOrientationType? cameraOrientation_)
        {
            SetFieldValue(3, 0, cameraOrientation_, Fit.SubfieldIndexMainField);
        }
        
        #endregion // Methods
    } // Class
} // namespace
