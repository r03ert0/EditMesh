//
//  MyDocument.h
//  EditMesh
//
//  Created by rOBERTO tORO on 13/11/2005.
//  Copyright __MyCompanyName__ 2005 . All rights reserved.
//


#import <Cocoa/Cocoa.h>
#import "MyMeshView.h"
#import "MyTextView.h"

@interface MyDocument : NSDocument
{
	IBOutlet MyMeshView			*view;
	IBOutlet MyTextView			*text;
	IBOutlet NSObjectController	*settings;
	MeshRec						M;
}
-(IBAction)importVertexData:(id)sender;
-(IBAction)importVertexParameterization:(id)sender;
-(IBAction)importTextureImage:(id)sender;
-(IBAction)exportVertexData:(id)sender;

-(IBAction)emChangeCentre:(id)sender;
-(IBAction)emDepth:(id)sender;
-(IBAction)emCurvature:(id)sender;			//
-(void)emIntegratedCurvature:(int)iter;
-(IBAction)emFlipTriangles:(id)sender;		//
-(IBAction)emParam2d:(id)sender;			//
-(IBAction)emRostrocaudal:(id)sender;		//
-(IBAction)emSmooth:(id)sender;				//
-(IBAction)emChangeCMapMinMax:(id)sender;
-(IBAction)emAddMesh:(id)sender;

-(void)addMesh:(char*)path;
@end
