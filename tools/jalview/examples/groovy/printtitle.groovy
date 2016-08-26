/*
 * Jalview - A Sequence Alignment Editor and Viewer (Version 2.9.0b2)
 * Copyright (C) 2015 The Jalview Authors
 * 
 * This file is part of Jalview.
 * 
 * Jalview is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *  
 * Jalview is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Jalview.  If not, see <http://www.gnu.org/licenses/>.
 * The Jalview Authors are detailed in the 'AUTHORS' file.
 */
// do something groovy in jalview
print "Hello World.\n";
def alf = Jalview.getAlignFrames();
for (ala in alf)
{
	// ala is an jalview.gui.AlignFrame object 
	print ala.getTitle()+"\n";
	// get the parent jalview.datamodel.Alignment from the alignment viewport
	def alignment = ala.viewport.alignment;
	// get the first sequence from the jalview.datamodel.Alignment object
	def seq = alignment.getSequenceAt(0); 
}
