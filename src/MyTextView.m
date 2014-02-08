#import "MyTextView.h"

@implementation MyTextView
- (id)initWithFrame:(NSRect)frameRect
{
	if ((self = [super initWithFrame:frameRect]) != nil)
	{
		app=nil;
	}
	return self;
}
-(void)keyDown:(NSEvent*)e
{
	//printf("key code:%i\n",[e keyCode]);
	if([e keyCode]==76||[e keyCode]==36)	// return and enter
	{
		if([e isARepeat]==NO)
		{
			unsigned int start,end;
			NSString	*s;
			[[self string] getLineStart:&start end:&end contentsEnd:nil forRange:(NSRange){[[self string] length]-1,1}];
			start+=2;
			s=[[self string] substringWithRange:(NSRange){start,end-start}];
			[history addObject:s];
			ihistory=0;
			[app applyScript:s];
		}
		
		NSRange range = { [[self string] length], 0 };
		[self setSelectedRange:range];
		[self insertText:@"\r> "];
	}
	else
	if([e keyCode]==126)	// up arrow
	{
		unsigned int	start,end;
		NSRange			range;
		
		ihistory++;
		if(ihistory>[history count])
			ihistory=[history count];
	
		[[self string] getLineStart:&start end:&end contentsEnd:nil forRange:(NSRange){[[self string] length]-1,1}];
		start+=2;
		range=(NSRange){start,end-start};
		[self setSelectedRange:range];
		[self insertText:[history objectAtIndex:[history count]-ihistory]];
	}
	else
	if([e keyCode]==125)	// down arrow
	{
		unsigned int	start,end;
		NSRange			range;
	
		[[self string] getLineStart:&start end:&end contentsEnd:nil forRange:(NSRange){[[self string] length]-1,1}];
		start+=2;
		range=(NSRange){start,end-start};
		[self setSelectedRange:range];
		
		ihistory--;
		if(ihistory<0)
		{
			ihistory=0;
			[self insertText:@""];
		}
		else
			[self insertText:[history objectAtIndex:[history count]-ihistory]];
	}
	else
	if([e keyCode]==48)	// tab
	{
		int				i,j;
		unsigned int	start,end;
		NSString		*s;
		NSRange			range;
		char			*first,*second;
		
		[[self string] getLineStart:&start end:&end contentsEnd:nil forRange:(NSRange){[[self string] length]-1,1}];
		start+=2;
		range=(NSRange){start,end-start};
		s=[[self string] substringWithRange:range];
		for(i=0;i<[history count];i++)
		{
			first=(char*)[s UTF8String];
			second=(char*)[[history objectAtIndex:[history count]-1-i] UTF8String];
			j=0;
			while(first[j]==second[j])
				j++;
			if(j==end-start)
			{
				[self setSelectedRange:range];
				[self insertText:[history objectAtIndex:[history count]-1-i]];
				break;
			}
		}
	}
	else
		[super keyDown:e];
}
#pragma mark -
-(void)setApp:(id)theApp
{
	app=theApp;
	history=[[NSMutableArray new] retain];
	ihistory=0;
	[self insertText:@"> "];
}
-(void)saveText:(NSString*)path
{
	[[self string] writeToFile:path atomically:YES];
}
@end
