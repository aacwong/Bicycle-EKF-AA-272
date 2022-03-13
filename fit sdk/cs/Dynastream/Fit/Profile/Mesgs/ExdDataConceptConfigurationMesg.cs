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
    /// Implements the ExdDataConceptConfiguration profile message.
    /// </summary>
    public class ExdDataConceptConfigurationMesg : Mesg
    {
        #region Fields
        #endregion

        /// <summary>
        /// Field Numbers for <see cref="ExdDataConceptConfigurationMesg"/>
        /// </summary>
        public sealed class FieldDefNum
        {
            public const byte ScreenIndex = 0;
            public const byte ConceptField = 1;
            public const byte FieldId = 2;
            public const byte ConceptIndex = 3;
            public const byte DataPage = 4;
            public const byte ConceptKey = 5;
            public const byte Scaling = 6;
            public const byte DataUnits = 8;
            public const byte Qualifier = 9;
            public const byte Descriptor = 10;
            public const byte IsSigned = 11;
            public const byte Invalid = Fit.FieldNumInvalid;
        }

        #region Constructors
        public ExdDataConceptConfigurationMesg() : base(Profile.GetMesg(MesgNum.ExdDataConceptConfiguration))
        {
        }

        public ExdDataConceptConfigurationMesg(Mesg mesg) : base(mesg)
        {
        }
        #endregion // Constructors

        #region Methods
        ///<summary>
        /// Retrieves the ScreenIndex field</summary>
        /// <returns>Returns nullable byte representing the ScreenIndex field</returns>
        public byte? GetScreenIndex()
        {
            Object val = GetFieldValue(0, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set ScreenIndex field</summary>
        /// <param name="screenIndex_">Nullable field value to be set</param>
        public void SetScreenIndex(byte? screenIndex_)
        {
            SetFieldValue(0, 0, screenIndex_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the ConceptField field</summary>
        /// <returns>Returns nullable byte representing the ConceptField field</returns>
        public byte? GetConceptField()
        {
            Object val = GetFieldValue(1, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set ConceptField field</summary>
        /// <param name="conceptField_">Nullable field value to be set</param>
        public void SetConceptField(byte? conceptField_)
        {
            SetFieldValue(1, 0, conceptField_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the FieldId field</summary>
        /// <returns>Returns nullable byte representing the FieldId field</returns>
        public byte? GetFieldId()
        {
            Object val = GetFieldValue(2, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set FieldId field</summary>
        /// <param name="fieldId_">Nullable field value to be set</param>
        public void SetFieldId(byte? fieldId_)
        {
            SetFieldValue(2, 0, fieldId_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the ConceptIndex field</summary>
        /// <returns>Returns nullable byte representing the ConceptIndex field</returns>
        public byte? GetConceptIndex()
        {
            Object val = GetFieldValue(3, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set ConceptIndex field</summary>
        /// <param name="conceptIndex_">Nullable field value to be set</param>
        public void SetConceptIndex(byte? conceptIndex_)
        {
            SetFieldValue(3, 0, conceptIndex_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the DataPage field</summary>
        /// <returns>Returns nullable byte representing the DataPage field</returns>
        public byte? GetDataPage()
        {
            Object val = GetFieldValue(4, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set DataPage field</summary>
        /// <param name="dataPage_">Nullable field value to be set</param>
        public void SetDataPage(byte? dataPage_)
        {
            SetFieldValue(4, 0, dataPage_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the ConceptKey field</summary>
        /// <returns>Returns nullable byte representing the ConceptKey field</returns>
        public byte? GetConceptKey()
        {
            Object val = GetFieldValue(5, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set ConceptKey field</summary>
        /// <param name="conceptKey_">Nullable field value to be set</param>
        public void SetConceptKey(byte? conceptKey_)
        {
            SetFieldValue(5, 0, conceptKey_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Scaling field</summary>
        /// <returns>Returns nullable byte representing the Scaling field</returns>
        public byte? GetScaling()
        {
            Object val = GetFieldValue(6, 0, Fit.SubfieldIndexMainField);
            if(val == null)
            {
                return null;
            }

            return (Convert.ToByte(val));
            
        }

        /// <summary>
        /// Set Scaling field</summary>
        /// <param name="scaling_">Nullable field value to be set</param>
        public void SetScaling(byte? scaling_)
        {
            SetFieldValue(6, 0, scaling_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the DataUnits field</summary>
        /// <returns>Returns nullable ExdDataUnits enum representing the DataUnits field</returns>
        public ExdDataUnits? GetDataUnits()
        {
            object obj = GetFieldValue(8, 0, Fit.SubfieldIndexMainField);
            ExdDataUnits? value = obj == null ? (ExdDataUnits?)null : (ExdDataUnits)obj;
            return value;
        }

        /// <summary>
        /// Set DataUnits field</summary>
        /// <param name="dataUnits_">Nullable field value to be set</param>
        public void SetDataUnits(ExdDataUnits? dataUnits_)
        {
            SetFieldValue(8, 0, dataUnits_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Qualifier field</summary>
        /// <returns>Returns nullable ExdQualifiers enum representing the Qualifier field</returns>
        public ExdQualifiers? GetQualifier()
        {
            object obj = GetFieldValue(9, 0, Fit.SubfieldIndexMainField);
            ExdQualifiers? value = obj == null ? (ExdQualifiers?)null : (ExdQualifiers)obj;
            return value;
        }

        /// <summary>
        /// Set Qualifier field</summary>
        /// <param name="qualifier_">Nullable field value to be set</param>
        public void SetQualifier(ExdQualifiers? qualifier_)
        {
            SetFieldValue(9, 0, qualifier_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the Descriptor field</summary>
        /// <returns>Returns nullable ExdDescriptors enum representing the Descriptor field</returns>
        public ExdDescriptors? GetDescriptor()
        {
            object obj = GetFieldValue(10, 0, Fit.SubfieldIndexMainField);
            ExdDescriptors? value = obj == null ? (ExdDescriptors?)null : (ExdDescriptors)obj;
            return value;
        }

        /// <summary>
        /// Set Descriptor field</summary>
        /// <param name="descriptor_">Nullable field value to be set</param>
        public void SetDescriptor(ExdDescriptors? descriptor_)
        {
            SetFieldValue(10, 0, descriptor_, Fit.SubfieldIndexMainField);
        }
        
        ///<summary>
        /// Retrieves the IsSigned field</summary>
        /// <returns>Returns nullable Bool enum representing the IsSigned field</returns>
        public Bool? GetIsSigned()
        {
            object obj = GetFieldValue(11, 0, Fit.SubfieldIndexMainField);
            Bool? value = obj == null ? (Bool?)null : (Bool)obj;
            return value;
        }

        /// <summary>
        /// Set IsSigned field</summary>
        /// <param name="isSigned_">Nullable field value to be set</param>
        public void SetIsSigned(Bool? isSigned_)
        {
            SetFieldValue(11, 0, isSigned_, Fit.SubfieldIndexMainField);
        }
        
        #endregion // Methods
    } // Class
} // namespace
