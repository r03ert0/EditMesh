/* MyTextView */

#import <Cocoa/Cocoa.h>

@interface MyTextView : NSTextView
{
	id app;
	NSMutableArray	*history;
	int				ihistory;
}
-(void)setApp:(id)theApp;
-(void)saveText:(NSString*)path;
@end
